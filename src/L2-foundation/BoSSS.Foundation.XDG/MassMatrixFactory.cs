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
using System.Runtime.CompilerServices;

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
        internal MassMatrixFactory(XDGSpaceMetrics __XDGSpaceMetrics) {
            XDGSpaceMetrics = __XDGSpaceMetrics;
            //this.MaxBasis = new Basis(XDGSpaceMetrics.GridDat, XDGSpaceMetrics.CutCellQuadOrder / 2); // bad choice;
            this.MaxBasis = new Basis(XDGSpaceMetrics.GridDat, 1);
            this.MaxTraceBasis = null;
        }



        /// <summary>
        /// the maximal supported basis if this factory
        /// </summary>
        public Basis MaxBasis {
            get;
            private set;
        }

        /// <summary>
        /// the maximal supported basis if this factory
        /// </summary>
        public Basis MaxTraceBasis {
            get;
            private set;
        }



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

        

        /// <summary>
        /// Provides access to a 'raw' form of the mass matrix, where only blocks for the cut cells
        /// are stored
        /// </summary>
        /// <returns>
        /// don't mess with those values
        /// </returns>
        public MassMatrixBlockContainer GetMassMatrixBlocks(Basis b, SpeciesId spc) {
            UpdateBlocks(b.Degree, new SpeciesId[] { spc });
            if (!b.IsSubBasis(this.MaxBasis))
                throw new NotSupportedException("requested basis exceeds maximally supported basis.");

            var MMB = this.MassBlocks[spc];
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

        /// <summary>
        /// cached mass-matrix blocks for cut-cells for TraceDG
        /// </summary>
        MassMatrixFactory.MassMatrixBlockContainer TraceMassBlocks;

        /// <summary>
        /// computes the mass matrices for a given mapping and accumulates the mass matrix to some other matrix.
        /// </summary>
        public void AccMassMatrix<T>(T M, UnsetteledCoordinateMapping mapping, IDictionary<SpeciesId, IEnumerable<double>> _alpha, bool inverse = false)
            where T : IMutableMatrixEx //
        {
            using (var tr = new FuncTrace()) {
                var _basisS = mapping.BasisS.ToArray();
                var ctx = _basisS[0].GridDat;
                int J = ctx.iLogicalCells.NoOfLocalUpdatedCells;
                if(!object.ReferenceEquals(this.XDGSpaceMetrics.GridDat, ctx))
                    throw new ArgumentException("grid object mismatch");
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

                
                tr.Info("Requesting Mass Matrix for degrees " + mapping.BasisS.Select(b => b.Degree).ToConcatString("", ", ", ";"));


                // compute the Mass-Blocks for the cut cells...
                // --------------------------------------------
                int _MaxDeg_vol = _basisS.Where(basis => !(basis is TraceDGBasis)).Select(basis => basis.Degree).Cat(-1).Max();
                int _MaxDeg_trc = _basisS.Where(basis => (basis is TraceDGBasis)).Select(basis => basis.Degree).Cat(-1).Max();
                var RequestedSpecies = _alpha.Keys;
                if(_MaxDeg_vol >= 0)
                    UpdateBlocks(_MaxDeg_vol, RequestedSpecies);
                if(_MaxDeg_trc >= 0)
                    UpdateTraceBlocks(_MaxDeg_trc);


                Dictionary<SpeciesId, MassMatrixBlockContainer>[] invBlocks = null;
                MassMatrixBlockContainer[] invTraceBlocks = null;
                if (inverse) {
                    if(_MaxDeg_vol >= 0)
                        invBlocks = GetInverseMassMatrixBlocks(RequestedSpecies, mapping.BasisS.Select(b => (b is XDGBasis) ? b.Degree : -1).ToArray());
                    if(_MaxDeg_trc >= 0)
                        invTraceBlocks = GetInverseTraceMassMatrixBlocks(mapping.BasisS.Select(b => (b is TraceDGBasis) ? b.Degree : -1).ToArray());
                }

                // set the matrix
                // --------------

                for (int iFld = 0; iFld < _basisS.Length; iFld++) { // loop over Variables in mapping...

                    // properties of the basis
                    XDGBasis Xbasis = _basisS[iFld] as XDGBasis;
                    TraceDGBasis traceBasis = _basisS[iFld] as TraceDGBasis;
                    Basis nonXbasis = Xbasis?.NonX_Basis ?? traceBasis?.NonX_Basis ?? _basisS[iFld];
                    Basis basis = _basisS[iFld];
                    int N = nonXbasis.Length;
                    if(Xbasis != null || object.ReferenceEquals(nonXbasis, _basisS[iFld])) {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // XDG basis OR standard DG basis
                        // (both have a volume-based mass matrix)
                        // Normally, the stadard DG basis has a mass matrix that is the identity matrix, i.e., the mass matrix for both phases sum up to be the identity matrix.
                        // There are, however, still codes like CNS-IBM which are programmed against older states of the XDG API, 
                        // they use a standard DG basis with a XDG configuration and therefore they (might) require this execution path.
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++



                        var regions = Xbasis?.Tracker?.Regions; // compute species indices with respect to **actual** (most recent timestep) regions !!!!
                                                                // species indices might have changed from one timestep to the next.
                                                         

                        if(traceBasis != null) {
                            throw new ApplicationException();
                        }
                        
                        foreach(var species in _alpha.Keys) { // loop over species...

                            double alpha = _alpha[species].ElementAt(iFld);
                            int[] jSub2jCell = this.MassBlocks[species].jSub2jCell;

                            if(alpha != 0.0) {

                                // get Mass-Matrix subblock
                                MultidimensionalArray Mass;
                                if(!inverse) {
                                    Mass = this.MassBlocks[species].MassMatrixBlocks;
                                } else {
                                    Mass = invBlocks[iFld][species].MassMatrixBlocks;
                                    Debug.Assert(Mass.GetLength(1) == N);
                                    Debug.Assert(Mass.GetLength(2) == N);
                                }

                                // set diagonal element 1.0 in all cells of the subgrid
                                // (later, we subtract 1.0 from the diagonal in cut cells)
                                var speciesMask = XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask(species);
                                {
                                    MultidimensionalArray eye = MultidimensionalArray.Create(N, N);
                                    eye.AccEye(1.0);

                                    foreach(int jCell in speciesMask.ItemEnum) {

                                        int iSpc = 0;
                                        if(Xbasis != null) {
                                            //int NoOfSpc = regions.GetNoOfSpecies(jCell);
                                            iSpc = regions.GetSpeciesIndex(species, jCell);
                                        }

                                        long n0 = mapping.GlobalUniqueCoordinateIndex(iFld, jCell, N * iSpc);
                                        M.AccBlock(n0, n0, alpha, eye); // might be faster than accessing M[n0 + n, n0 + n] sequentially
                                        //for(int n = 0; n < N; n++) {
                                        //    M[n0 + n, n0 + n] += alpha;
                                        //}

                                    }
                                }

                                // now, set then Non-Identity blocks in the cut cells
                                var speciesBitMask = speciesMask.GetBitMask();
                                {
                                    int JSUB = jSub2jCell.Length;

                                    for(int jsub = 0; jsub < JSUB; jsub++) { // loop over cut cells

                                        int jCell = jSub2jCell[jsub];
                                        if(!speciesBitMask[jCell])
                                            continue;


                                        int iSpc = 0;
                                        if(Xbasis != null) {
                                            //int NoOfSpc = regions.GetNoOfSpecies(jCell);
                                            iSpc = regions.GetSpeciesIndex(species, jCell);
                                        }

                                        var MassSub = Mass.ExtractSubArrayShallow(new int[] { jsub, 0, 0 }, new int[] { jsub - 1, N - 1, N - 1 });

                                        long i0 = mapping.GlobalUniqueCoordinateIndex(iFld, jCell, N * iSpc);

                                        //var M2 = M.CloneAs();
                                        //var MassSub2 = MassSub.CloneAs();

                                        var DiagBkup = MassSub.GetDiagonal();
                                        MassSub.AccEye(-1.0); // subtract the diagonal that we added previously
                                        M.AccBlock(i0, i0, alpha, MassSub);
                                        MassSub.SetDiagonal(DiagBkup); // restore diagonal (alternatively, .AccEye(+alpha) could be used but that might introduce round-off errors?);

                                        //for(int n = 0; n < N; n++) { // loop over rows (within block)
                                        //    for(int m = 0; m < N; m++) { // loop over columns (within block)
                                        //        M[i0 + n, i0 + m] += alpha * MassSub[n, m] - (m == n ? alpha : 0.0);
                                        //    }
                                        //}

                                        //var dist = M2.MatrixDistFrobenius(M);
                                        //var distb = MassSub2.FrobeniusDistance(MassSub);
                                        //Console.WriteLine(jsub + "of" + JSUB + ": " +  dist + "   " + distb + "   " + alpha);
                                    }
                                }
                            }
                        }
                    } else if(traceBasis != null) {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // Trace DG basis 
                        // (bases on a surface integral, zero outside of cut cells)
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        nonXbasis = traceBasis.NonX_Basis;
                       

                        double alpha = _alpha.First().Value.ElementAt(iFld);
                        int[] jSub2jCell = this.TraceMassBlocks.jSub2jCell;

                        if(alpha != 0.0) {

                            // get Mass-Matrix subblock
                            MultidimensionalArray Mass;
                            if(!inverse) {
                                Mass = this.TraceMassBlocks.MassMatrixBlocks;
                            } else {
                                Mass = invTraceBlocks[iFld].MassMatrixBlocks;
                                Debug.Assert(Mass.GetLength(1) == N);
                                Debug.Assert(Mass.GetLength(2) == N);
                            }

                            // now, set then Non-Identity blocks in the cut cells
                            //var speciesBitMask = speciesMask.GetBitMask();
                            {
                                int JSUB = jSub2jCell.Length;

                                for(int jsub = 0; jsub < JSUB; jsub++) { // loop over cut cells

                                    int jCell = jSub2jCell[jsub];

                                    var MassSub = Mass.ExtractSubArrayShallow(new int[] { jsub, 0, 0 }, new int[] { jsub - 1, N - 1, N - 1 });
                                    long i0 = mapping.GlobalUniqueCoordinateIndex(iFld, jCell, 0);
                                    M.AccBlock(i0, i0, alpha, MassSub);

                                    
                                }
                            }
                        }

                    } else {
                        throw new NotImplementedException($"unkonwn basis class {_basisS[iFld]}");
                    }
                    


                }
            }
        }

        private void UpdateBlocks(int _MaxDeg, IEnumerable<SpeciesId> RequestedSpecies) {
            using (var tr = new FuncTrace()) {
                if(_MaxDeg < 0)
                    return;

                if (MassBlocks == null)
                    MassBlocks = new Dictionary<SpeciesId, MassMatrixFactory.MassMatrixBlockContainer>();

                // Try to update the blocks only once: for the Basis with highest Polynomial Degree
                if (_MaxDeg > this.MaxBasis.Degree) {
                    tr.Info("Mass Matrix requested for degree: " + _MaxDeg);
                    MassBlocks.Clear();
                    this.MaxBasis = new Basis(this.MaxBasis.GridDat, _MaxDeg);
                } else {
                    tr.Info("Mass Matrix for basis of degree: " + this.MaxBasis.Degree);
                }
                Basis nonXbasis = this.MaxBasis;

                // compute Blocks
                {
                    var SpeciesToDo = RequestedSpecies.Except(MassBlocks.Keys);
                    if (SpeciesToDo.Count() > 0) {
                        // only for species that we haven't done/cached yet
                        Dictionary<SpeciesId, MassMatrixFactory.MassMatrixBlockContainer> _MassBlocks;
                        ComputeMassMatrixBlocks(SpeciesToDo, out _MassBlocks, nonXbasis, this.XDGSpaceMetrics);
                        this.MassBlocks.AddRange(_MassBlocks);
                    }
                }
            }

        }

        private void UpdateTraceBlocks(int _MaxDeg) {
            using(var tr = new FuncTrace()) {
                if(_MaxDeg < 0)
                    return;


                // // Try to update the blocks only once: for the Basis with highest Polynomial Degree
                if(this.MaxTraceBasis == null || _MaxDeg > this.MaxTraceBasis.Degree) {
                    tr.Info("Mass Matrix requested for degree: " + _MaxDeg);
                    TraceMassBlocks = null;
                    this.MaxTraceBasis = new Basis(this.MaxBasis.GridDat, _MaxDeg);
                } else {
                    tr.Info("Mass Matrix for basis of degree: " + this.MaxBasis.Degree);
                }
                Basis nonXbasis = this.MaxTraceBasis;

                // compute Blocks
                ComputeTraceMassMatrixBlocks(out this.TraceMassBlocks, nonXbasis, this.XDGSpaceMetrics);
            }

        }



        /// <summary>
        /// 
        /// </summary>
        public class MassMatrixBlockContainer {

            /// <summary>
            /// Mass matrix blocks
            /// - 1st index: subgrid cell index (which map to local cell indices by <see cref="jSub2jCell"/>)
            /// - 2nd index: diagonal block row;
            /// - 3rd index: diagonal block column;
            /// </summary>
            public MultidimensionalArray MassMatrixBlocks;
       
            /// <summary>
            /// mapping from sub-grid index (correlates to 1st index of <see cref="MassMatrixBlocks"/>) to local cell index
            /// </summary>
            public int[] jSub2jCell;


            /// <summary>
            /// mapping from cell index to subgrid index (correlates to 1st index of <see cref="MassMatrixBlocks"/>)
            /// </summary>
            public int[] jCell2jSub;

            /// <summary>
            /// The integration domain for mass matrices; this contains only cut cells.
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
                tracer.Info("Mass Matrix order: " + b.Degree + " -> dim = " + b.Length);

                int quadorder = homie.CutCellQuadOrder;


                // define domains and allocate memory
                // ==================================


                foreach (var Species in _SpeciesIds) { // loop over species...

                    // interation dom
                    var _IntegrationDomain = homie.LevelSetRegions.GetSpeciesMask(Species).Intersect(homie.LevelSetRegions.GetCutCellMask());

                    // domain for mass-matrix blocks 
                    var _BlockDomain = _IntegrationDomain;

                    // alloc mem for blocks
                    var _MassMatrixBlocksSpc = MultidimensionalArray.Create(_BlockDomain.NoOfItemsLocally, Nnx, Nnx);

                    // Subgrid index to cell index
                    int[] _jSub2jCell = _BlockDomain.ItemEnum.ToArray();

                    // the inverse:
                    int[] _jCell2jSub = new int[ctx.iLogicalCells.Count];
                    _jCell2jSub.SetAll(int.MinValue);
                    for(int i = 0; i < _jSub2jCell.Length; i++) {
                        _jCell2jSub[_jSub2jCell[i]] = i;
                    }

                    // cell to subgrid index
                    Result.Add(Species, new MassMatrixBlockContainer() {
                        IntegrationDomain = _IntegrationDomain,
                        MassMatrixBlocks = _MassMatrixBlocksSpc,
                        jCell2jSub = _jCell2jSub,
                        jSub2jCell = _jSub2jCell
                    });

                }

                // compute blocks
                // ==============

                foreach (var Species in _SpeciesIds) {
                    int[] jSub2jCell = Result[Species].jSub2jCell;
                    int[] jCell2jSub = Result[Species].jCell2jSub;
                    // get quad scheme
                    CellQuadratureScheme scheme = schemeHelper.GetVolumeQuadScheme(Species, IntegrationDomain: Result[Species].IntegrationDomain);
                    
                    // result storage
                    var MassMatrixBlocksSpc = Result[Species].MassMatrixBlocks;
                    //for(int i = 0; i < jSub2jCell.Length; i++) {
                    //    for(int n = 0; n < Nnx; n++)
                    //        MassMatrixBlocksSpc[i,n,n] = 0.0;
                    //}

                    tracer.Info("mass matrix quad order: " + quadorder);

                    // compute the products of the basis functions:
                    //int BlockCnt = -1;
                    CellMask speciesCells = homie.LevelSetRegions.GetSpeciesMask(Species);
                    var compRule = scheme.Compile(ctx, quadorder);
                    var quad = CellQuadrature.GetQuadrature(
                        new int[] { Nnx, Nnx },
                        ctx,
                        compRule,
                        delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                            // Del_Evaluate
                            // ~~~~~~~~~~~~~
                            /*
                            if(QR.Nodes.NoOfNodes > 10000) {
                                var phi = homie.Tracker.LevelSets[0] as SinglePhaseField;
                                Console.WriteLine("Rule with: " + QR.NoOfNodes + " nodes");
                                for(int j = i0; j < i0+Length; j++) {
                                    Console.WriteLine("in cell: " + j);
                                    Console.WriteLine(homie.Tracker.GridDat.Grid.Cells[j].ToString());
                                    homie.Tracker.GridDat.Grid.Cells[j].TransformationParams.SaveToStream(Console.Out);
                                    Console.WriteLine("Degree of Level-Set: " + phi.Basis.Degree);
                                    Console.WriteLine("DG coordinates: " + phi.Coordinates.GetRow(j).ToConcatString("", ", ", ";"));
                                }
                            }
                            */

                            var BasisVal = b.CellEval(QR.Nodes, i0, Length);
                            EvalResult.Multiply(1.0, BasisVal, BasisVal, 0.0, "ikmn", "ikm", "ikn");

                        },
                        delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                            // Del_SaveIntegrationResults
                            // ~~~~~~~~~~~~~~~~~~~~~~~~~~

                            for(int i = 0; i < Length; i++) {
                                int jCell = i0 + i;
                                //BlockCnt++;
                                int BlockIdx = jCell2jSub[jCell];

                                // insert ID block in agglomeration target cells (if necessary):
                                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                var Block = MassMatrixBlocksSpc.ExtractSubArrayShallow(BlockIdx, -1, -1);
                                //while(BlockCell[BlockCnt] < jCell) {
                                //    // agglomeration source/target cell that is not cut
                                //    // mass matrix is identity (full) or zero (void)
                                //    Block.Clear();
                                //    if(speciesCells.Contains(BlockCell[BlockCnt])) {
                                //        // cell is full
                                //        for(int nn = 0; nn < Nnx; nn++) {
                                //            Block[nn, nn] = 1.0;
                                //        }
                                //    }
                                //    BlockCnt++;
                                //    Block = MassMatrixBlocksSpc.ExtractSubArrayShallow(BlockCnt, -1, -1);
                                //}

                                // store computed block
                                // - - - - - - - - - - - 
                                Debug.Assert(jSub2jCell[BlockIdx] == jCell);
                                MassMatrixBlocksSpc.ExtractSubArrayShallow(BlockIdx, -1, -1)
                                    .Set(ResultsOfIntegration.ExtractSubArrayShallow(i, -1, -1));
#if DEBUG
                                for (int n = 0; n < Nnx; n++)
                                    for (int m = 0; m < Nnx; m++)
                                        Debug.Assert(Block[n, m] == Block[m, n]);
#endif
                            }
                        });
                    quad.Execute();
                    // ------------------------------------ quadrature end.
                }
            }
        }

        static internal void ComputeTraceMassMatrixBlocks(
            out MassMatrixBlockContainer Result,
            Basis b,
            XDGSpaceMetrics homie) {

            using(var tracer = new FuncTrace()) {
                if(b is XDGBasis)
                    throw new ArgumentException();
                var ctx = homie.GridDat;

                var schemeHelper = homie.XQuadSchemeHelper;
                int Nnx = b.Length;
                tracer.Info("Mass Matrix order: " + b.Degree + " -> dim = " + b.Length);

                int quadorder = homie.CutCellQuadOrder;


                // define domains and allocate memory
                // ==================================


                { 
                    
                    // interation dom
                    var _IntegrationDomain = homie.LevelSetRegions.GetCutCellMask();

                    // domain for mass-matrix blocks 
                    var _BlockDomain = _IntegrationDomain;

                    // alloc mem for blocks
                    var _MassMatrixBlocksSpc = MultidimensionalArray.Create(_BlockDomain.NoOfItemsLocally, Nnx, Nnx);

                    // Subgrid index to cell index
                    int[] _jSub2jCell = _BlockDomain.ItemEnum.ToArray();

                    // the inverse:
                    int[] _jCell2jSub = new int[ctx.iLogicalCells.Count];
                    _jCell2jSub.SetAll(int.MinValue);
                    for(int i = 0; i < _jSub2jCell.Length; i++) {
                        _jCell2jSub[_jSub2jCell[i]] = i;
                    }

                    // cell to subgrid index
                    Result = new MassMatrixBlockContainer() {
                        IntegrationDomain = _IntegrationDomain,
                        MassMatrixBlocks = _MassMatrixBlocksSpc,
                        jCell2jSub = _jCell2jSub,
                        jSub2jCell = _jSub2jCell
                    };

                }

                // compute blocks
                // ==============

                {
                    int[] jSub2jCell = Result.jSub2jCell;
                    int[] jCell2jSub = Result.jCell2jSub;
                    // get quad scheme
                    CellQuadratureScheme scheme = schemeHelper.GetLevelSetquadScheme(0, homie.SpeciesList.First(), Result.IntegrationDomain);

                    // result storage
                    var MassMatrixBlocksSpc = Result.MassMatrixBlocks;
                    tracer.Info("trace mass matrix quad order: " + quadorder);

                    // compute the products of the basis functions:
                    var compRule = scheme.Compile(ctx, quadorder);
                    var quad = CellQuadrature.GetQuadrature(
                        new int[] { Nnx, Nnx },
                        ctx,
                        compRule,
                        delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                            // Del_Evaluate
                            // ~~~~~~~~~~~~~
                            var BasisVal = b.CellEval(QR.Nodes, i0, Length);
                            EvalResult.Multiply(1.0, BasisVal, BasisVal, 0.0, "ikmn", "ikm", "ikn");

                        },
                        delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                            // Del_SaveIntegrationResults
                            // ~~~~~~~~~~~~~~~~~~~~~~~~~~

                            for(int i = 0; i < Length; i++) {
                                int jCell = i0 + i;
                                //BlockCnt++;
                                int BlockIdx = jCell2jSub[jCell];

                                // insert ID block in agglomeration target cells (if necessary):
                                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                var Block = MassMatrixBlocksSpc.ExtractSubArrayShallow(BlockIdx, -1, -1);
                               
                                // store computed block
                                // - - - - - - - - - - - 
                                Debug.Assert(jSub2jCell[BlockIdx] == jCell);
                                MassMatrixBlocksSpc.ExtractSubArrayShallow(BlockIdx, -1, -1)
                                    .Set(ResultsOfIntegration.ExtractSubArrayShallow(i, -1, -1));
#if DEBUG
                                for(int n = 0; n < Nnx; n++)
                                    for(int m = 0; m < Nnx; m++)
                                        Debug.Assert(Block[n, m] == Block[m, n]);
#endif
                            }
                        });
                    quad.Execute();
                    // ------------------------------------ quadrature end.
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

        MassMatrixFactory.MassMatrixBlockContainer[] GetInverseTraceMassMatrixBlocks(int[] Degrees) {
            int[] Ns = Degrees.Select(p => this.MaxBasis.Polynomials[0].Where(poly => poly.AbsoluteDegree <= p).Count()).ToArray();
            for(int iKref = 1; iKref < this.MaxBasis.Polynomials.Count; iKref++) {
                int[] _Ns = Degrees.Select(p => this.MaxBasis.Polynomials[0].Where(poly => poly.AbsoluteDegree <= p).Count()).ToArray();
                if(!ArrayTools.ListEquals(Ns, _Ns))
                    throw new NotSupportedException();
            }

            MassMatrixBlockContainer[] InvMassBlocks = new MassMatrixBlockContainer[Degrees.Length];


            for(int i = 0; i < Degrees.Length; i++) {
                int pSame = -1;
                //pSame = Array.IndexOf(Degrees, Degrees[i], 0, i);
                for(int j = 0; j < i; j++) {
                    if((Degrees[j] == Degrees[i])
                        //&& (VariableAgglomerationSwitch[j] == VariableAgglomerationSwitch[i])
                        ) {
                        pSame = j;
                        break;
                    }
                }

                if(pSame >= 0) {
                    InvMassBlocks[i] = InvMassBlocks[pSame];
                } else {
                    InvMassBlocks[i] = InvertMassMatrixBlock(Ns[i], this.TraceMassBlocks);
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


                    MassMatrixBlocksInv.Add(Species,
                        InvertMassMatrixBlock(N, MassMatrixBlocks[Species]));
                }
            }

        }

        private static MassMatrixBlockContainer InvertMassMatrixBlock(int N, MassMatrixBlockContainer BlockContainer) {
            var mblk = BlockContainer.MassMatrixBlocks;
            //var mblkB4agg = MassMatrixBlocks[Species].MassMatrixBlocks_B4Agglom;
            var invmblk = MultidimensionalArray.Create(mblk.GetLength(0), N, N);

            var ret = new MassMatrixBlockContainer() {
                MassMatrixBlocks = invmblk,
                //jCell2jSub = MassMatrixBlocks[Species].jCell2jSub,
                jSub2jCell = BlockContainer.jSub2jCell,
                IntegrationDomain = BlockContainer.IntegrationDomain
            };


            int B = mblk.GetLength(0);
            Debug.Assert(mblk.GetLength(2) >= N);

            MultidimensionalArray blk = MultidimensionalArray.Create(N, N);
            MultidimensionalArray inv_Blk = MultidimensionalArray.Create(N, N);


            // loop over all cut cells (mass matrix blocks)...
            for(int b = 0; b < B; b++) {
                MultidimensionalArray _blk;
                //if ((aggSw == true) || (mblkB4agg[b] == null)) {
                _blk = mblk.ExtractSubArrayShallow(new int[] { b, 0, 0 }, new int[] { b - 1, N - 1, N - 1 });
                //} else {
                //    Debug.Assert(aggSw == false);
                //    Debug.Assert(mblkB4agg[b] != null);
                //    _blk = mblkB4agg[b].ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { N - 1, N - 1 });
                //}

                blk.Set(_blk);
                if(blk.ContainsNanOrInf(true, true))
                    throw new ArithmeticException("NAN/INF in mass matrix.");

                if(blk.InfNorm() == 0.0) {
                    // a zero-block
                    invmblk.ExtractSubArrayShallow(b, -1, -1).Clear();
                } else {
                    //blk.Invert(inv_Blk);
                    try {
                        inv_Blk.Clear();
                        inv_Blk.Acc(1.0, blk);
                        inv_Blk.InvertSymmetrical();

#if DEBUG
                        for(int n = 0; n < N; n++)
                            for(int m = 0; m < N; m++)
                                Debug.Assert(inv_Blk[n, m] == inv_Blk[m, n]);

#endif
                    } catch(ArithmeticException) {
                        //Console.WriteLine("WARNING: indefinite mass matrix.");
                        blk.InvertTo(inv_Blk);
#if DEBUG
                        for(int n = 0; n < N; n++)
                            for(int m = 0; m < N; m++)
                                Debug.Assert(Math.Abs(inv_Blk[n, m] - inv_Blk[m, n]) / (Math.Abs(inv_Blk[n, m]) + Math.Abs(inv_Blk[m, n])) <= 1.0e-9);
#endif
                    }

                    if(blk.ContainsNanOrInf(true, true))
                        throw new ArithmeticException("NAN/INF in inverse mass matrix.");

                    invmblk.SetSubMatrix(inv_Blk, b, -1, -1);
                }
            }


            return ret;
        }
    }
}
