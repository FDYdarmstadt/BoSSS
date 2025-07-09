using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using ilPSP.Tracing;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Microsoft.CodeAnalysis.CSharp.Syntax;
using System.Collections;
using BoSSS.Foundation.Comm;
using System.Diagnostics;
using ilPSP.Utils;

namespace BoSSS.Solution.AdvancedSolvers {
    /// <summary>
    /// Aggregation basis for trace fields on level-set surfaces
    /// </summary>
    public class TraceDGAggregationBasis : AggregationGridBasis {

        /// <summary>
        /// Trace DG basis on original grid
        /// </summary>
        public TraceDGBasis TraceDGBasis {
            get;
            private set;
        }

        private int[] m_NoOfSpecies;
        
        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="tb">
        /// Trace DG basis on original grid
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
        internal TraceDGAggregationBasis(TraceDGBasis tb, TraceDGAggregationBasis parentBasis, AggregationGridData ag, MultidimensionalArray[] inj)
            : base(tb.NonX_Basis, parentBasis, ag, inj) {
            using(new FuncTrace()) {
                this.TraceDGBasis = tb;
            }
        }

        /// <summary>
        /// Updates species mapping for trace fields
        /// </summary>
        /// <param name="Agglomerator"></param>
        public void Update(MultiphaseCellAgglomerator Agglomerator) {
            using(new FuncTrace()) {
                UpdateSpeciesMapping(Agglomerator);
                m_XCompositeBasis = null;
            }
        }


        private void UpdateBaseGridInjector(int jagg) {
            var LsTrk = this.TraceDGBasis.Tracker;
            int N = this.DGBasis.Length;
            int[][] compCells = this.AggGrid.iLogicalCells.AggregateCellToParts;



            // mass matrix in the aggregate cell, before we orthonormalize:
            var AggCellMMb4Ortho = MultidimensionalArray.Create(N, N);

            int NoOfSpc_jagg = this.m_NoOfSpecies[jagg];
            int[] AggCell = compCells[jagg];
            int K = AggCell.Length;

            {
               
                { 
                    bool partiallyEmptyCell = false; // true: at least one of the base cells which form 
                                                     //       aggregate cell 'jagg' is empty with respect to the TraceDG-basis
                                                     //         => this will alter the projection operator.
                    bool Empty = true; 
                    for(int k = 0; k < K; k++) {
                        if(TraceDGBasis.GetLength(AggCell[k]) <= 0) {
                            partiallyEmptyCell = true;
                            Empty = false;
                        }
                    }

                    if(Empty) {
                        this.m_XCompositeBasis[jagg] = MultidimensionalArray.Create(K, N, N);

                    } else if(partiallyEmptyCell) {

                        // compute mass matrix in aggregate cell 'jagg' for species index 'iSpc_agg'
                        // -------------------------------------------------------------------------

                        AggCellMMb4Ortho.Clear();
                        for(int k = 0; k < K; k++) { // loop over the base cells in the aggregate cell
                            if(TraceDGBasis.GetLength(AggCell[k]) > 0) {

                                var ExPolMtx = base.GetCompositeBasis(jagg).ExtractSubArrayShallow(k, -1, -1);
                                AggCellMMb4Ortho.Multiply(1.0, ExPolMtx, ExPolMtx, 1.0, "lm", "im", "il");
                            }
                        }


                        // change to orthonormal basis
                        // ---------------------------
                        MultidimensionalArray B = MultidimensionalArray.Create(N, N);
                        try {
                            AggCellMMb4Ortho.SymmetricLDLInversion(B, default(double[]));
                        } catch(ArithmeticException ae) {
                            #region diagnostic_output                                     
                            //Console.Error.WriteLine("ArithmeticException in XdgAggregationBasis.Update() at MG level " + this.AggGrid.MgLevel + " for aggregate cell " + jagg + " and species index " + iSpc_agg);
                            //continue;  
                            Console.Error.WriteLine(ae.GetType() + ": " + ae.Message);
                            Console.Error.WriteLine("Mesh level: " + AggGrid.MgLevel);
                            Console.Error.WriteLine("Aggregate cell " + jagg);
                            int[] parts = (AggGrid.iLogicalCells?.AggregateCellToParts[jagg]) ?? new int[0];
                            Console.Error.WriteLine("Aggregate Cell: " + parts.ToConcatString("{", ",", "}"));
                            var LevSet0 = LsTrk.LevelSets[0] as LevelSet;
                            var Marker = new SinglePhaseField(new Basis(LevSet0.GridDat, 0), "marker");
                            BitArray ba = new BitArray(Marker.GridDat.iLogicalCells.NoOfLocalUpdatedCells);
                            foreach(var j in parts) {
                                Marker.SetMeanValue(j, 1.0);
                                ba[j] = true;
                            }

                            

                            for(int k = 0; k < K; k++) { // loop over the base cells in the aggregate cell
                                {

                                    var ExPolMtx = base.GetCompositeBasis(jagg).ExtractSubArrayShallow(k, -1, -1);
                                    ExPolMtx.SaveToTextFile("Expol-" + k + ".txt");
                                    var Mama = MultidimensionalArray.Create(ExPolMtx.NoOfRows, ExPolMtx.NoOfCols);
                                    Mama.Multiply(1.0, ExPolMtx, ExPolMtx, 1.0, "lm", "im", "il");



                                    bool IsPosDef;
                                    try {
                                        Mama.Cholesky();
                                        IsPosDef = true;
                                    } catch(ArithmeticException) {
                                        IsPosDef = false;
                                    }

                                    Console.Error.WriteLine($"Part {k}: posdef? {IsPosDef}");
                                }
                            }

                            AggCellMMb4Ortho.SaveToTextFile("indef.txt");
                            throw ae;
                            #endregion
                        }

                        this.m_XCompositeBasis[jagg] = MultidimensionalArray.Create(K, N, N);

                        var X_ExPolMtx = this.m_XCompositeBasis[jagg];
                        var NonX_ExPolMtx = GetCompositeBasis(jagg);
                        X_ExPolMtx.Allocate(NonX_ExPolMtx.Lengths); // should not reallocate if lengths stay the same;

                        X_ExPolMtx.Multiply(1.0, NonX_ExPolMtx, B, 0.0, "imn", "imk", "kn");
                    } else {
                        // no empty cell with respect to species 'iSpc_agg' in aggregate cell
                        // --> non-XDG projector can be used.
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        this.m_XCompositeBasis[jagg] = null;
                    }
                }
            }
        }


        const bool PerformAgglomeration = false;


        private void UpdateSpeciesMapping(MultiphaseCellAgglomerator Agglomerator) {
            using(new FuncTrace()) {
                var LsTrk = this.TraceDGBasis.Tracker;
                var agg = base.AggGrid;
                int[][] compCells = agg.iLogicalCells.AggregateCellToParts;
                int JAGG = compCells.Length;
                Debug.Assert(JAGG == agg.iLogicalCells.Count);
                int JAGGup = agg.iLogicalCells.NoOfLocalUpdatedCells;
                
                // Get cut cell mask from level set tracker
                var CutCellMask = LsTrk.Regions.GetCutCellMask().GetBitMask();

                if(m_NoOfSpecies == null) {
                    m_NoOfSpecies = new int[compCells.Length];
                } else {
                    Debug.Assert(m_NoOfSpecies.Length == JAGG);
                }

                // for TraceDG, we agglomerate if ANY species is agglomerated (shall be ALL species, might differ from the actually used ones)
                BitArray agglomeratedCells = null;
                if(PerformAgglomeration) {
#pragma warning disable 162
                    foreach(var SpId in Agglomerator.SpeciesList) {
                        var _agglomeratedCellsSpecies = Agglomerator.GetAgglomerator(SpId).AggInfo.SourceCells.GetBitMaskWithExternal();
                        if(agglomeratedCells == null)
                            agglomeratedCells = _agglomeratedCellsSpecies;
                        else
                            agglomeratedCells = agglomeratedCells.Or(_agglomeratedCellsSpecies);
                    }
#pragma warning restore 162
                }

                // loop over all aggregate (multigrid) cells...
                for(int jagg = 0; jagg < JAGGup; jagg++) {
                    int[] compCell = compCells[jagg];
                    int K = compCell.Length;

                    // check no of species:
                    bool hasCut = false;
                    bool isAgglom = PerformAgglomeration;
                    foreach(int j in compCell) {
                        hasCut = hasCut || CutCellMask[j];
                        isAgglom = isAgglom && agglomeratedCells[j];
                        if(hasCut == true && isAgglom == false)
                            break; // perf. opt; no need to look at the rest.
                    }


                    m_NoOfSpecies[jagg] = (hasCut && !isAgglom) ? 1 : 0; // contains a species if any of the base cells contains a cut cell and is not agglomerated
                } 
                
                // MPI exchange for parallel consistency
                m_NoOfSpecies.MPIExchange(agg);
            }
        }


        MultidimensionalArray[] m_XCompositeBasis;
        BitArray m_XCompositeBasisUpdated;


        public MultidimensionalArray GetXCompositeBasis(int jAgg, int iSpcAgg) {
            if(m_XCompositeBasis == null) {
                m_XCompositeBasis = new MultidimensionalArray[base.AggGrid.iLogicalCells.NoOfLocalUpdatedCells];
                m_XCompositeBasisUpdated = new BitArray(m_XCompositeBasis.Length);
            }

            if(!m_XCompositeBasisUpdated[jAgg]) {
                UpdateBaseGridInjector(jAgg);
                m_XCompositeBasisUpdated[jAgg] = true;
            }

            MultidimensionalArray Trf;
            if(this.m_XCompositeBasis[jAgg] == null)
                Trf = base.GetCompositeBasis(jAgg);
            else
                Trf = this.m_XCompositeBasis[jAgg];

            return Trf;
        }




        /// <summary>
        /// Number of species in composite/aggregate cell <paramref name="jAgg"/>.
        /// Only species with non-zero measure are counted.
        /// </summary>
        override public int GetNoOfSpecies(int jAgg) {
            return m_NoOfSpecies[jAgg];
        }



        /// <summary>
        /// Maximum length of basis
        /// </summary>
        public override int MaximalLength {
            get {
                return this.TraceDGBasis.MaximalLength;
            }
        }

        /// <summary>
        /// Minimal length of basis
        /// </summary>
        public override int MinimalLength {
            get {
                return this.TraceDGBasis.MinimalLength;
            }
        }

        /// <summary>
        /// Gets length for given cell and polynomial degree
        /// </summary>
        public override int GetLength(int jCell, int p) {
            return m_NoOfSpecies[jCell] * base.GetLength(jCell, p);
        }

        /// <summary>
        /// Gets maximal length for polynomial degree; equal to a standard <see cref="AggregationGridBasis"/>, since 
        /// a TraceDG-field only has a single value in cut-cells (in contrast to XDG fields, which can have multiple values)
        /// </summary>
        public override int GetMaximalLength(int p) {
            return base.GetMaximalLength(p);
        }

        /// <summary>
        /// Gets minimal length for polynomial degree; 0, since TraceDG-fields are undefined outside of cut-cells
        /// </summary>
        public override int GetMinimalLength(int p) {
            return 0;
        }
    }
}
