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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Creation of XDG mass-matrices.
    /// </summary>
    public class MassMatrixFactory {

        /// <summary>
        /// owner object.
        /// </summary>
        public XDGSpaceMetrics XDGSpaceMetrics {
            get;
            private set;
        }


        /// <summary>
        /// ctor.
        /// </summary>
        public MassMatrixFactory(XDGSpaceMetrics __XDGSpaceMetrics) {
            XDGSpaceMetrics = __XDGSpaceMetrics;
            this.MaxBasis = new Basis(XDGSpaceMetrics.GridDat, XDGSpaceMetrics.CutCellQuadOrder / 2);
        }

        /*
        /// <summary>
        /// quadrature order used for creating the volume rule
        /// </summary>
        public int QuadOrder {
            get {
                return m_quadorder;
            }
        }
        */

        /// <summary>
        /// the maximal supported basis if this factory
        /// </summary>
        public Basis MaxBasis {
            get;
            private set;
        }

        //Basis m_MaxNonXBasis {
        //    get {
        //        if (MaxBasis is XDGBasis) {
        //            return ((XDGBasis)MaxBasis).NonX_Basis;
        //        } else {
        //            return MaxBasis;
        //        }
        //    }
        //}


        //MultiphaseCellAgglomerator m_agglomerator;

        //LevelSetTracker m_LsTrk;

   

        //public XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant {
        //    get;
        //    private set;
        //}

        ///// <summary>
        ///// Computes the mass matrices for a given mapping, for all species in <see cref="AvailableSpecies"/>
        ///// </summary>
        //public BlockDiagonalMatrix GetMassMatrix_depr(UnsetteledCoordinateMapping mapping, bool inverse) {
        //    double[] alpha = new double[mapping.BasisS.Count];
        //    alpha.SetAll(1.0);
        //    return GetMassMatrix_depr(mapping, alpha, inverse, this.AvailableSpecies.ToArray());
        //}

        ///// <summary>
        ///// Computes the mass matrices for a given mapping.
        ///// </summary>
        //public BlockDiagonalMatrix GetMassMatrix_depr(UnsetteledCoordinateMapping mapping, double[] _alpha, bool inverse, params SpeciesId[] Spc) {
        //    Dictionary<SpeciesId, IEnumerable<double>> alpha = new Dictionary<SpeciesId, IEnumerable<double>>();
        //    foreach (var species in Spc)
        //        alpha.Add(species, _alpha.CloneAs());
        //    return GetMassMatrix_depr(mapping, alpha, inverse);
        //}

        /// <summary>
        /// computes the mass matrices for a given mapping.
        /// </summary>
        /// <param name="mapping">
        /// </param>
        /// <param name="alpha">
        /// 'Fine-grained' scaling for blocks, for each species and each variable in the <paramref name="mapping"/>.
        /// The default value of null maps to 1.0 for each species and variable.
        /// </param>
        /// <param name="inverse">
        /// Return the inverse mass matrix.
        /// </param>
        /// <param name="agg">
        /// Required if an agglomerated Mass matrix is requested.
        /// </param>
        public BlockMsrMatrix GetMassMatrix(UnsetteledCoordinateMapping mapping, IDictionary<SpeciesId, IEnumerable<double>> alpha = null, bool inverse = false) {
            if (alpha == null) {
                int NoVar = mapping.NoOfVariables;
                alpha = new Dictionary<SpeciesId, IEnumerable<double>>();
                double[] AllOne = new double[NoVar];
                AllOne.SetAll(1.0);
                foreach (var spc in this.XDGSpaceMetrics.SpeciesList) {
                    alpha.Add(spc, AllOne);
                }
            }
            var Return = new BlockMsrMatrix(mapping, mapping);
            AccMassMatrix(Return, mapping, alpha, inverse);
            return Return;
        }

        /// <summary>
        /// Computes the mass matrices for a given mapping, for all species in <see cref="AvailableSpecies"/>
        /// </summary>
        public BlockMsrMatrix GetMassMatrix(UnsetteledCoordinateMapping mapping, bool inverse) {
            double[] alpha = new double[mapping.BasisS.Count];
            alpha.SetAll(1.0);
            return GetMassMatrix(mapping, alpha, inverse, this.AvailableSpecies.ToArray());
        }

        /// <summary>
        /// Computes the mass matrices for a given mapping.
        /// </summary>
        public BlockMsrMatrix GetMassMatrix(UnsetteledCoordinateMapping mapping, double[] _alpha, bool inverse, params SpeciesId[] Spc) {
            Dictionary<SpeciesId, IEnumerable<double>> alpha = new Dictionary<SpeciesId, IEnumerable<double>>();
            foreach (var species in Spc)
                alpha.Add(species, _alpha.CloneAs());
            return GetMassMatrix(mapping, alpha, inverse);
        }

        /*
        /// <summary>
        /// computes the mass matrices for a given mapping.
        /// </summary>
        /// <param name="mapping">
        /// </param>
        /// <param name="alpha">
        /// 'Fine-grained' scaling for blocks, for each species and each variable in the <paramref name="mapping"/>.
        /// The default value of null maps to 1.0 for each species and variable.
        /// </param>
        /// <param name="VariableAgglomerationSwitch">
        /// Switch to turn agglomeration on/off for each variable; default=null=on.
        /// </param>
        /// <param name="inverse">
        /// Return the inverse mass matrix.
        /// </param>
        public BlockDiagonalMatrix GetMassMatrix_depr(UnsetteledCoordinateMapping mapping, IDictionary<SpeciesId, IEnumerable<double>> alpha = null, bool inverse = false, bool[] VariableAgglomerationSwitch = null) {
            if (alpha == null) {
                int NoVar = mapping.NoOfVariables;
                alpha = new Dictionary<SpeciesId, IEnumerable<double>>();
                double[] AllOne = new double[NoVar];
                AllOne.SetAll(1.0);
                foreach (var spc in this.m_LsTrk.SpeciesIdS) {
                    alpha.Add(spc, AllOne);
                }
            }
            var Return = new BlockDiagonalMatrix(mapping.LocalLength, mapping.MaxTotalNoOfCoordinatesPerCell);
            AccMassMatrix(Return, mapping, alpha, inverse, VariableAgglomerationSwitch);
            return Return;
        }
        */

        /// <summary>
        /// Provides access to a 'raw' form of the mass matrix, where only blocks for the cut cells
        /// are stored
        /// </summary>
        /// <returns>
        /// don't mess with those values
        /// </returns>
        public MassMatrixBlockContainer GetMassMatrixBlocks(Basis b, SpeciesId spc) {
            if (!b.IsSubBasis(this.MaxBasis))
                throw new NotSupportedException("requested basis exceeds maximally supported basis.");

            UpdateBlocks(b.Degree, new SpeciesId[] { spc });

            var MMB = this.MassBlocks[spc];
//#if DEBUG
//            {
//                var CCBit = XDGSpaceMetrics.LevelSetRegions.GetCutCellMask().GetBitMask();
//                //var AggCells = this.Agglomerator.GetAgglomerator(spc).AggInfo.SourceCells;
//                //var AggCellsBit = AggCells.GetBitMask();
//                //var AggAffectedBit = this.Agglomerator.GetAgglomerator(spc).AggInfo.AllAffectedCells.GetBitMaskWithExternal();

//                int J = XDGSpaceMetrics.GridDat.Cells.NoOfLocalUpdatedCells;

//                int JSUB = MMB.jSub2jCell.Length;
//                for (int jSub = 0; jSub < JSUB; jSub++) {
//                    var Block = MMB.MassMatrixBlocks.ExtractSubArrayShallow(jSub, -1, -1);
//                    int jCell = MMB.jSub2jCell[jSub];
//                    double BlockNorm = Block.InfNorm();

//                    //if (jCell < J) {
//                    //    Debug.Assert(CCBit[jCell] || AggAffectedBit[jCell]);
//                    //    Debug.Assert(AggCellsBit[jCell] == (BlockNorm == 0));
//                    //}
//                }

//                //foreach (var jCellAgg in AggCells.ItemEnum) {
//                //    int jSub = MMB.jCell2jSub[jCellAgg];
//                //    var Block = MMB.MassMatrixBlocks.ExtractSubArrayShallow(jSub, -1, -1);
//                //    double BlockNorm = Block.InfNorm();

//                //    Debug.Assert(BlockNorm == 0.0);
//                //}
//            }

//#endif
            return MMB;
        }

        /// <summary>
        /// a collection of all species that are available 
        /// </summary>
        public IEnumerable<SpeciesId> AvailableSpecies {
            get {
                return this.XDGSpaceMetrics.SpeciesList;
            }
        }

        /// <summary>
        /// cached mass-matrix blocks for cut-cells
        /// </summary>
        Dictionary<SpeciesId, MassMatrixFactory.MassMatrixBlockContainer> MassBlocks;

        ///// <summary>
        ///// cached inverse mass-matrix blocks for cut-cells
        ///// </summary>
        //Dictionary<SpeciesId, MassMatrixFactory.MassMatrixBlockContainer> InverseMassBlocks;

        

        /// <summary>
        /// computes the mass matrices for a given mapping and accumulates the mass matrix to some other matrix.
        /// </summary>
        public void AccMassMatrix<T>(T M, UnsetteledCoordinateMapping mapping, IDictionary<SpeciesId, IEnumerable<double>> _alpha, bool inverse = false)
            where T : IMutableMatrixEx //
        {
            using (new FuncTrace()) {
                var _basisS = mapping.BasisS.ToArray();
                var ctx = _basisS[0].GridDat;
                int J = ctx.iLogicalCells.NoOfLocalUpdatedCells;

                //if (VariableAgglomerationSwitch == null) {
                //    VariableAgglomerationSwitch = new bool[mapping.BasisS.Count];
                //    VariableAgglomerationSwitch.SetAll(true);
                //} else {
                //    if (VariableAgglomerationSwitch.Length != mapping.BasisS.Count)
                //        throw new ArgumentException();
                //}

                //if(VariableAgglomerationSwitch.Any() && (agg == null)) {
                //    throw new ArgumentException("Cell Agglomerator is required.");
                //}

                if (!M.RowPartitioning.EqualsPartition(mapping))
                    throw new ArgumentException("Mismatch in row mapping.");
                if (!M.ColPartition.EqualsPartition(mapping))
                    throw new ArgumentException("Mismatch in column mapping.");

                if (!_alpha.Keys.IsSubsetOf(this.AvailableSpecies))
                    throw new ArgumentException("trying to obtain mass matrix for species that is not supported by this factory.");

                foreach (var __alpha in _alpha.Values) {
                    if (__alpha.Count() != _basisS.Length)
                        throw new ArgumentException("Number of alpha's must match number of variables/fields in mapping.");
                }

                LevelSetTracker.LevelSetRegions regions = null;// XDGSpaceMetrics.LevelSetRegions;


                // compute the Mass-Blocks for the cut cells...
                // --------------------------------------------
                int _MaxDeg = _basisS.Max(b => b.Degree);
                var RequestedSpecies = _alpha.Keys;
                UpdateBlocks(_MaxDeg, RequestedSpecies);

                Dictionary<SpeciesId, MassMatrixBlockContainer>[] invBlocks = null;
                if (inverse) {
                    invBlocks = GetInverseMassMatrixBlocks(RequestedSpecies, mapping.BasisS.Select(b => b.Degree).ToArray());
                    //VariableAgglomerationSwitch = new bool[VariableAgglomerationSwitch.Length];
                    //VariableAgglomerationSwitch.SetAll(true); // we already took care about this in the 'GetInverseMassMatrixBlocks'
                }

                // set the matrix
                // --------------

                for (int fld = 0; fld < _basisS.Length; fld++) { // loop over Variables in mapping...

                    // loop over species...
                    foreach (var species in _alpha.Keys) {

                        double alpha = _alpha[species].ElementAt(fld);
                        int[] jSub2jCell = this.MassBlocks[species].jSub2jCell;

                        if (alpha != 0.0) {

                            // properties of the basis
                            XDGBasis Xbasis = _basisS[fld] as XDGBasis;
                            Basis nonXbasis;
                            Basis basis = _basisS[fld];
                            if (Xbasis != null) {
                                nonXbasis = Xbasis.NonX_Basis;
                                regions = Xbasis.Tracker.Regions; // compute species indices with respect to **actual** regions !!!!
                                //if (!object.ReferenceEquals(Xbasis.Tracker, this.XDGSpaceMetrics.))
                                //    throw new ArgumentException();
                            } else {
                                nonXbasis = _basisS[fld];
                            }
                            int N = nonXbasis.Length;

                            // get Mass-Matrix subblock
                            MultidimensionalArray Mass;
                            if (!inverse) {
                                Mass = this.MassBlocks[species].MassMatrixBlocks;
                            } else {
                                Mass = invBlocks[fld][species].MassMatrixBlocks;
                                Debug.Assert(Mass.GetLength(1) == N);
                                Debug.Assert(Mass.GetLength(2) == N);
                            }

                            // set diagonal element 1.0 in all cells of the subgrid
                            // (later, we subtract 1.0 from the diagonal in cut cells)
                            var speciesMask = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(species);
                            {

                                foreach (int jCell in speciesMask.ItemEnum) {

                                    int iSpc = 0;
                                    if (Xbasis != null) {
                                        int NoOfSpc = regions.GetNoOfSpecies(jCell);
                                        iSpc = regions.GetSpeciesIndex(species, jCell);
                                    }

                                    int n0 = mapping.GlobalUniqueCoordinateIndex(fld, jCell, N * iSpc);

                                    for (int n = 0; n < N; n++) {
                                        M[n0 + n, n0 + n] += alpha;
                                    }

                                }
                            }

                            // now, set then Non-Identity blocks in the cut cells
                            var speciesBitMask = speciesMask.GetBitMask();
                            {
                                int JSUB = jSub2jCell.Length;

                                for (int jsub = 0; jsub < JSUB; jsub++) { // loop over cut cells

                                    int jCell = jSub2jCell[jsub];
                                    if (!speciesBitMask[jCell])
                                        continue;


                                    int iSpc = 0;
                                    if (Xbasis != null) {
                                        int NoOfSpc = regions.GetNoOfSpecies(jCell);
                                        iSpc = regions.GetSpeciesIndex(species, jCell);
                                    }

                                    var MassSub = Mass.ExtractSubArrayShallow(new int[] { jsub, 0, 0 }, new int[] { jsub - 1, N - 1, N - 1 });
                                    //MultidimensionalArray MassSub;
                                    //if (VariableAgglomerationSwitch[fld] || (Mass_B4Agglom[jsub] == null)) {
                                    //    // block with agglomeration
                                    //    MassSub = Mass.ExtractSubArrayShallow(new int[] { jsub, 0, 0 }, new int[] { jsub - 1, N - 1, N - 1 });
                                    //} else {
                                    //    // block without agglomeration

                                    //    Debug.Assert(VariableAgglomerationSwitch[fld] == false);
                                    //    Debug.Assert(Mass_B4Agglom[jsub] != null);
                                    //    MassSub = Mass_B4Agglom[jsub].ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { N - 1, N - 1 });
                                    //}


                                    int i0 = mapping.GlobalUniqueCoordinateIndex(fld, jCell, N * iSpc);

                                    for (int n = 0; n < N; n++) { // loop over rows (within block)
                                        for (int m = 0; m < N; m++) { // loop over columns (within block)
                                            M[i0 + n, i0 + m] += alpha * MassSub[n, m] - (m == n ? alpha : 0.0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                //// cell-aglomeration
                //// -----------------

                //if (inverse == true)
                //    //throw new ApplicationException("todo: double-check with agglomeration");
                //    Console.WriteLine("todo: double-check with agglomeration");

                //this.m_agglomerator.ManipulateMassMatrix_Mk2(M, mapping);
            }
        }

        private void UpdateBlocks(int _MaxDeg, IEnumerable<SpeciesId> RequestedSpecies) {
            if (MassBlocks == null)
                MassBlocks = new Dictionary<SpeciesId, MassMatrixFactory.MassMatrixBlockContainer>();

            // ..., but only once: for the Basis with highest Polynomial Degree
            if (_MaxDeg > this.MaxBasis.Degree)
                throw new ArgumentException();

            Basis nonXbasis = this.MaxBasis;

            // compute Blocks
            {
                var SpeciesToDo = RequestedSpecies.Except(MassBlocks.Keys);
                if (SpeciesToDo.Count() > 0) {
                    // only for species that we haven't done/cached yet
                    Dictionary<SpeciesId, MassMatrixFactory.MassMatrixBlockContainer> _MassBlocks;
                    ComputeMassMatrixBlocks(SpeciesToDo, out _MassBlocks, nonXbasis, this.XDGSpaceMetrics); // m_quadorder, m_LsTrk, this.m_agglomerator, MomentFittingVariant);
                    this.MassBlocks.AddRange(_MassBlocks);
                }
            }


        }

        /// <summary>
        /// 
        /// </summary>
        /// <remarks>
        /// Note that due to cell agglomeration,
        /// there may be non-identity mass matrices in uncut cells:
        /// this happens when a cut cell is agglomerated to an un-cut cell.
        /// </remarks>
        public class MassMatrixBlockContainer {

            /// <summary>
            /// Mass matrix blocks, with agglomeration applied. 
            /// 1st index: subgrid cell index (which map to local cell indices by <see cref="jSub2jCell"/>)
            /// 2nd index: diagonal block row;
            /// 3rd index: diagonal block column;
            /// </summary>
            public MultidimensionalArray MassMatrixBlocks;

            /*
            /// <summary>
            /// 'Backup' of mass matrix blocks before agglomeration; if an entry is null, the corresponding cell is not affected by agglomeration.
            /// </summary>
            public MultidimensionalArray[] MassMatrixBlocks_B4Agglom;
            */

            /// <summary>
            /// mapping from sub-grid index (correlates to 1st index of <see cref="MassMatrixBlocks"/>) to local cell index
            /// </summary>
            public int[] jSub2jCell;

            ///// <summary>
            ///// the inverse mapping of <see cref="jSub2jCell"/>
            ///// </summary>
            //public Dictionary<int, int> jCell2jSub;


            /// <summary>
            /// The integration domain for mass matrices; this contains only cut cells.
            /// It may not correspond to <see cref="jSub2jCell"/>, since 
            /// some cut-cells may be agglomerated to uncut cells.
            /// </summary>
            internal CellMask IntegrationDomain;

        }

        static internal void ComputeMassMatrixBlocks(
            IEnumerable<SpeciesId> _SpeciesIds,
            out Dictionary<SpeciesId, MassMatrixBlockContainer> Result,
            Basis b,
            XDGSpaceMetrics homie) {

            using (var tracer = new FuncTrace()) {
                if (b is XDGBasis)
                    throw new ArgumentException();
                var ctx = homie.GridDat;

                Result = new Dictionary<SpeciesId, MassMatrixBlockContainer>();
                var schemeHelper = homie.XQuadSchemeHelper;
                int Nnx = b.Length;

                int quadorder = homie.CutCellQuadOrder;


                // define domains and allocate memory
                // ==================================


                foreach (var Species in _SpeciesIds) { // loop over species...

                    // interation dom
                    var _IntegrationDomain = homie.LevelSetRegions.GetSpeciesMask(Species).Intersect(homie.LevelSetRegions.GetCutCellMask());

                    // domain for mass-matrix blocks (include agglomeration targets)
                    var _BlockDomain = _IntegrationDomain; //.Union(Agg.GetAgglomerator(Species).AggInfo.AllAffectedCells);

                    // alloc mem for blocks
                    var _MassMatrixBlocksSpc = MultidimensionalArray.Create(_BlockDomain.NoOfItemsLocally, Nnx, Nnx);

                    // Subgrid index to cell index
                    int[] _jSub2jCell = _BlockDomain.ItemEnum.ToArray();

                    // cell to subgrid index
                    //Dictionary<int, int> _jCell2jSub;
                    //if (Agg.GetAgglomerator(Species).AggInfo.AgglomerationPairs.Length > 0) {
                    //    _jCell2jSub = new Dictionary<int, int>();
                    //    for (int i = 0; i < _jSub2jCell.Length; i++) {
                    //        _jCell2jSub.Add(_jSub2jCell[i], i);
                    //    }
                    //} else {
                    //    _jCell2jSub = null;
                    //}

                    Result.Add(Species, new MassMatrixBlockContainer() {
                        IntegrationDomain = _IntegrationDomain,
                        MassMatrixBlocks = _MassMatrixBlocksSpc,
                        //jCell2jSub = _jCell2jSub,
                        jSub2jCell = _jSub2jCell
                    });

                }

                // compute blocks
                // ==============

                foreach (var Species in _SpeciesIds) {
                    // get quad scheme
                    CellQuadratureScheme scheme = schemeHelper.GetVolumeQuadScheme(Species, IntegrationDomain: Result[Species].IntegrationDomain);

                    // result storage
                    var MassMatrixBlocksSpc = Result[Species].MassMatrixBlocks;

                    tracer.Info("mass matrix quad order: " + quadorder);

                    // compute the products of the basis functions:
                    int BlockCnt = -1;
                    int[] BlockCell = Result[Species].jSub2jCell;
                    CellMask speciesCells = homie.LevelSetRegions.GetSpeciesMask(Species);
                    CellQuadrature.GetQuadrature(
                        new int[] { Nnx, Nnx },
                        ctx,
                        scheme.Compile(ctx, quadorder),
                        delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                            // Del_Evaluate
                            // ~~~~~~~~~~~~~
                            var BasisVal = b.CellEval(QR.Nodes, i0, Length);
                            EvalResult.Multiply(1.0, BasisVal, BasisVal, 0.0, "ikmn", "ikm", "ikn");

                        },
                        delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                            // Del_SaveIntegrationResults
                            // ~~~~~~~~~~~~~~~~~~~~~~~~~~

                            for (int i = 0; i < Length; i++) {
                                int jCell = i0 + i;
                                BlockCnt++;

                                // insert ID block in agglomeration target cells (if necessary):
                                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                var Block = MassMatrixBlocksSpc.ExtractSubArrayShallow(BlockCnt, -1, -1);
                                while (BlockCell[BlockCnt] < jCell) {
                                    // agglomeration source/target cell that is not cut
                                    // mass matrix is identity (full) or zero (void)
                                    Block.Clear();
                                    if (speciesCells.Contains(BlockCell[BlockCnt]))
                                    {
                                        // cell is full
                                        for (int nn = 0; nn < Nnx; nn++) {
                                            Block[nn, nn] = 1.0;
                                        }
                                    }
                                    BlockCnt++;
                                    Block = MassMatrixBlocksSpc.ExtractSubArrayShallow(BlockCnt, -1, -1);
                                }

                                // store computed block
                                // - - - - - - - - - - - 
                                Debug.Assert(BlockCell[BlockCnt] == jCell);
                                MassMatrixBlocksSpc.ExtractSubArrayShallow(BlockCnt, -1, -1)
                                    .Set(ResultsOfIntegration.ExtractSubArrayShallow(i, -1, -1));
#if DEBUG
                                for (int n = 0; n < Nnx; n++)
                                    for (int m = 0; m < Nnx; m++)
                                        Debug.Assert(Block[n, m] == Block[m, n]);
#endif
                            }
                        }).Execute();
                    // ------------------------------------ quadrature end.

                    BlockCnt++;
                    while (BlockCnt < MassMatrixBlocksSpc.GetLength(0)) {
                        // agglomeration source/target cell that is not cut
                        // mass matrix is identity (full) or zero (void)
                        var Block = MassMatrixBlocksSpc.ExtractSubArrayShallow(BlockCnt, -1, -1);
                        Block.Clear();
                        if (speciesCells.Contains(BlockCell[BlockCnt])) 
                        {
                            // cell is full
                            for (int nn = 0; nn < Nnx; nn++) {
                                Block[nn, nn] = 1.0;
                            }
                        }
                        BlockCnt++;
                    }


                    /*
                    // test mass matrix for positive definiteness
                    {
                        int JSUB = MassMatrixBlocksSpc.GetLength(0);
                        SubGrid Idom = null;

                        int failCount = 0;
                        var PosDefiniteTest = new FullMatrix(Nnx, Nnx);

                        for (int jsub = 0; jsub < JSUB; jsub++) {
                            PosDefiniteTest.Clear();
                            PosDefiniteTest.Acc(MassMatrixBlocksSpc.ExtractSubArrayShallow(jsub, -1, -1), 1.0);

                            try {
                                PosDefiniteTest.Clear();
                                PosDefiniteTest.Acc(MassMatrixBlocksSpc.ExtractSubArrayShallow(jsub, -1, -1), 1.0);
                                PosDefiniteTest.InvertSymmetrical();

                                //PosDefiniteTest.Clear();
                                //PosDefiniteTest.AccEye(1.0);
                                //PosDefiniteTest.Acc(MassMatrixBlocksSpc.ExtractSubArrayShallow(jsub, -1, -1), -1.0);
                                //PosDefiniteTest.InvertSymmetrical();
                            } catch (ArithmeticException ae) {
                                if (Idom == null)
                                    Idom = new SubGrid(scheme.Domain);

                                int jCell = Idom.SubgridIndex2LocalCellIndex[jsub];
                                long Gid = Tracker.GridDat.Cells.GetCell(jCell).GlobalID;


                                double volFrac = Tracker.GetSpeciesVolume(jCell, Species)/ctx.Cells.GetCellVolume(jCell);

                                var errString = string.Format("Indefinite mass matrix in cell: globalId = {0}, local index = {1}, species {2}; \n   cell volume fraction: {3};\n   [{4}]", Gid, jCell, Tracker.GetSpeciesName(Species), volFrac, ae.Message);
                                tracer.Logger.Error(errString);
                                //Console.WriteLine(errString);
                                failCount++;
                            }
                        }

                        if (failCount > 0) {
                            var errString = string.Format("Indefinite mass matrix in {0} of {1} cut cells", failCount, JSUB);
                            tracer.Logger.Error(errString);
                            Console.WriteLine(errString);
                        } else {
                            Console.WriteLine("No indefinite mass matrix blocks");
                        }

                    }
                    // */

                    // backup before agglomeration (required if we wanna treat e.g. velocity in DG and pressure in XDG)
                    //MultidimensionalArray[] massMatrixBlocksB4Agglom = new MultidimensionalArray[Result[Species].jSub2jCell.Length];
                    //Result[Species].MassMatrixBlocks_B4Agglom = massMatrixBlocksB4Agglom;
                    //var _jCell2jSub = Result[Species].jCell2jSub;
                    //int J = ctx.Cells.NoOfLocalUpdatedCells;
                    //foreach (var pair in Agg.GetAgglomerator(Species).AggInfo.AgglomerationPairs) {

                    //    foreach (int jCell in new int[] { pair.jCellSource, pair.jCellTarget }) { // create a backup of source and target cell
                    //        if (jCell >= J)
                    //            continue;

                    //        int jSub = _jCell2jSub[jCell];

                    //        if (massMatrixBlocksB4Agglom[jSub] == null) {
                    //            massMatrixBlocksB4Agglom[jSub] = MassMatrixBlocksSpc.ExtractSubArrayShallow(jSub, -1, -1).CloneAs();
                    //        }
                    //    }
                    //}

                    // agglomeration
                    //Agg.GetAgglomerator(Species).ManipulateMassMatrixBlocks(MassMatrixBlocksSpc, b, Result[Species].jSub2jCell, Result[Species].jCell2jSub);
                    //throw new NotImplementedException("todo");
                }
            }
        }

        Dictionary<SpeciesId, MassMatrixFactory.MassMatrixBlockContainer>[] GetInverseMassMatrixBlocks(IEnumerable<SpeciesId> RequestedSpecies, int[] Degrees) {
            int[] Ns = Degrees.Select(p => this.MaxBasis.Polynomials[0].Where(poly => poly.AbsoluteDegree <= p).Count()).ToArray();
            for(int iKref = 1; iKref < this.MaxBasis.Polynomials.Count; iKref++ ) {
                int[] _Ns = Degrees.Select(p => this.MaxBasis.Polynomials[0].Where(poly => poly.AbsoluteDegree <= p).Count()).ToArray();
                if(!ArrayTools.ListEquals(Ns, _Ns))
                    throw new NotSupportedException();
            }

            Dictionary<SpeciesId, MassMatrixFactory.MassMatrixBlockContainer>[] InvMassBlocks = new Dictionary<SpeciesId, MassMatrixBlockContainer>[Degrees.Length];

            for (int i = 0; i < Degrees.Length; i++) {
                int pSame = -1;
                //pSame = Array.IndexOf(Degrees, Degrees[i], 0, i);
                for (int j = 0; j < i; j++) {
                    if ((Degrees[j] == Degrees[i]) 
                        //&& (VariableAgglomerationSwitch[j] == VariableAgglomerationSwitch[i])
                        ) {
                        pSame = j;
                        break;
                    }
                }

                if (pSame >= 0) {
                    InvMassBlocks[i] = InvMassBlocks[pSame];
                } else {
                    InvertMassMatrixBlocks(out InvMassBlocks[i], this.MassBlocks, Ns[i]);
                }
            }

            return InvMassBlocks;
        }

        static internal void InvertMassMatrixBlocks(
            out Dictionary<SpeciesId, MassMatrixBlockContainer> MassMatrixBlocksInv,
            Dictionary<SpeciesId, MassMatrixBlockContainer> MassMatrixBlocks,
            int N) {
            using (new FuncTrace()) {

                MassMatrixBlocksInv = new Dictionary<SpeciesId, MassMatrixBlockContainer>();
                foreach (var Species in MassMatrixBlocks.Keys) {

                    var mblk = MassMatrixBlocks[Species].MassMatrixBlocks;
                    //var mblkB4agg = MassMatrixBlocks[Species].MassMatrixBlocks_B4Agglom;
                    var invmblk = MultidimensionalArray.Create(mblk.GetLength(0), N, N);
                    MassMatrixBlocksInv.Add(Species,
                        new MassMatrixBlockContainer() {
                            MassMatrixBlocks = invmblk,
                            //jCell2jSub = MassMatrixBlocks[Species].jCell2jSub,
                            jSub2jCell = MassMatrixBlocks[Species].jSub2jCell,
                            IntegrationDomain = MassMatrixBlocks[Species].IntegrationDomain
                        });

                    int B = mblk.GetLength(0);
                    Debug.Assert(mblk.GetLength(2) >= N);

                    MultidimensionalArray blk = MultidimensionalArray.Create(N, N);
                    MultidimensionalArray inv_Blk = MultidimensionalArray.Create(N, N);


                    // loop over all cut cells (mass matrix blocks)...
                    for (int b = 0; b < B; b++) {
                        MultidimensionalArray _blk;
                        //if ((aggSw == true) || (mblkB4agg[b] == null)) {
                        _blk = mblk.ExtractSubArrayShallow(new int[] { b, 0, 0 }, new int[] { b - 1, N - 1, N - 1 });
                        //} else {
                        //    Debug.Assert(aggSw == false);
                        //    Debug.Assert(mblkB4agg[b] != null);
                        //    _blk = mblkB4agg[b].ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { N - 1, N - 1 });
                        //}

                        blk.Set(_blk);
                        if (blk.ContainsNanOrInf(true, true))
                            throw new ArithmeticException("NAN/INF in mass matrix.");

                        if (blk.InfNorm() == 0.0) {
                            // a zero-block
                            invmblk.ExtractSubArrayShallow(b, -1, -1).Clear();
                        } else {
                            //blk.Invert(inv_Blk);
                            try {
                                inv_Blk.Clear();
                                inv_Blk.Acc(1.0, blk);
                                inv_Blk.InvertSymmetrical();

#if DEBUG
                                for (int n = 0; n < N; n++)
                                    for (int m = 0; m < N; m++)
                                        Debug.Assert(inv_Blk[n, m] == inv_Blk[m, n]);

#endif
                            } catch (ArithmeticException) {
                                Console.WriteLine("WARNING: indefinite mass matrix.");
                                blk.InvertTo(inv_Blk);
#if DEBUG
                                for (int n = 0; n < N; n++)
                                    for (int m = 0; m < N; m++)
                                        Debug.Assert(Math.Abs(inv_Blk[n, m] - inv_Blk[m, n]) / (Math.Abs(inv_Blk[n, m]) + Math.Abs(inv_Blk[m, n])) <= 1.0e-9);
#endif
                            }

                            if (blk.ContainsNanOrInf(true, true))
                                throw new ArithmeticException("NAN/INF in inverse mass matrix.");

                            invmblk.SetSubMatrix(inv_Blk, b, -1, -1);
                        }

                        /*
                        bool Choleski_failed = false;
                        bool Symmelim_failed = false;
                        try {
                            blk.Initialize(_blk);
                            blk.InvertSymmetrical();
                        } catch (ArithmeticException) {
                            Choleski_failed = true;
                        }
                        try {
                            blk.Initialize(_blk);
                            FullMatrix.TiredRoutine(blk);
                        } catch (ArithmeticException) {
                            Symmelim_failed = true;
                        }

                        if (Choleski_failed || Symmelim_failed) {
                            Console.WriteLine("Mass matrix defect ({0},{1}): species {2}, jsub = {3};", Choleski_failed, Symmelim_failed,  LsTrk.GetSpeciesName(Species), b);
                        
                            int J = LsTrk.Ctx.Grid.NoOfUpdateCells;
                            if (MassErrors == null) {
                                MassErrors = new Dictionary<SpeciesId, System.Collections.BitArray>();
                            }
                            if (!MassErrors.ContainsKey(Species)) {
                                MassErrors.Add(Species, new System.Collections.BitArray(J));
                            }
                            var _MassErrors = MassErrors[Species];

                            int[] globalIdx = cutted.SubgridIndex2LocalCellIndex;
                            int jCell = globalIdx[b];

                            _MassErrors[jCell] = true;

                        }
                        //*/
                    }


                }
            }

        }
    }
}
