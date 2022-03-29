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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// algebraic cell agglomeration to avoid small cut-cells;
    /// </summary>
    public class MultiphaseCellAgglomerator : ICutCellMetrics {

        /// <summary>
        /// number of agglomerations in all species
        /// </summary>
        public int TotalNumberOfAgglomerations {
            get;
            private set;
        }

        /// <summary>
        /// threshold factor that determines when cells are agglomerated
        /// </summary>
        public double AgglomerationThreshold {
            get;
            private set;
        }

        /// <summary>
        /// similar to <see cref="AgglomerationThreshold"/>, but for old timesteps
        /// </summary>
        public double[] AgglomerationThreshold_Oldtimesteps {
            get;
            private set;
        }


        /// <summary>
        /// The quadrature order used for computing cell volumes and edge areas.
        /// </summary>
        public int CutCellQuadratureOrder {
            get {
                return NonAgglomeratedMetrics.CutCellQuadratureOrder;
            }
        }

        /// <summary>
        /// The kind of HMF which is be used for computing cell volumes.
        /// </summary>
        public XQuadFactoryHelper.MomentFittingVariants HMFvariant {
            get {
                return NonAgglomeratedMetrics.HMFvariant;
            }
        }

        /// <summary>
        /// All species for which agglomeration is available.
        /// </summary>
        public IEnumerable<SpeciesId> SpeciesList {
            get {
                return NonAgglomeratedMetrics.SpeciesList;
            }
        }

        /// <summary>
        /// Link to the tracker.
        /// </summary>
        public LevelSetTracker Tracker {
            get;
            private set;
        }

        /// <summary>
        /// Cut-cell length scales before agglomeration
        /// </summary>
        public CutCellMetrics NonAgglomeratedMetrics {
            get;
            private set;
        }

        public XDGSpaceMetrics XDGSpaceMetrics {
            get;
            private set;
        }


        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="__AgglomerationTreshold">see <see cref="AgglomerationThreshold"/></param>
        /// <param name="oldTs__AgglomerationTreshold">
        /// Agglomeration thresholds for  _previous_ timesteps.
        ///  The number of entries in this array determines how many previous timesteps are considered
        /// (<see cref="LevelSetTracker.HistoryLength"/>).
        /// </param>
        /// <param name="AgglomerateNewborn">
        /// 'Newborn' cut-cells 
        /// are agglomerated to cells which already exist in the previous  timestep.
        /// </param>
        /// <param name="AgglomerateDecased">
        /// 'Deceased' cut-cells are agglomerated to cells which exist in the next timestep.
        /// </param>
        /// <param name="ExceptionOnFailedAgglomeration">
        /// If true, an exception is thrown for 
        /// any cell which should be agglomerated, if no neighbour is found.
        /// </param>
        /// <param name="CutCellsQuadOrder">
        /// cut-cell quadrature order for the quadrature rule that is used to determine cell volumes;
        /// </param>
        /// <param name="Spc"></param>
        /// <param name="lsTrk"></param>
        /// <param name="NewbornAndDecasedThreshold">
        /// Volume fraction threshold at which a cut-cell counts as newborn, resp. deceased, see <paramref name="AgglomerateNewborn"/>, <paramref name="AgglomerateDecased"/>;
        /// this should typically be the same order which is used to evaluate the XDG operator matrix.
        /// </param>
        internal MultiphaseCellAgglomerator(
            LevelSetTracker lsTrk,
            SpeciesId[] Spc, int CutCellsQuadOrder,
            double __AgglomerationTreshold,
            bool AgglomerateNewborn = false, bool AgglomerateDecased = false, bool ExceptionOnFailedAgglomeration = true,
            double[] oldTs__AgglomerationTreshold = null,
            double NewbornAndDecasedThreshold = 1.0e-6
            ) {
            MPICollectiveWatchDog.Watch();
            if (__AgglomerationTreshold < 0.0 || __AgglomerationTreshold >= 1.0)
                throw new ArgumentOutOfRangeException();

            if (NewbornAndDecasedThreshold < 0.0 || NewbornAndDecasedThreshold >= 1.0)
                throw new ArgumentOutOfRangeException();

            this.Tracker = lsTrk;

            this.XDGSpaceMetrics = lsTrk.GetXDGSpaceMetrics(Spc, CutCellsQuadOrder, 1);
            this.NonAgglomeratedMetrics = lsTrk.GetXDGSpaceMetrics(Spc, CutCellsQuadOrder, 1).CutCellMetrics;
            this.AgglomerationThreshold = __AgglomerationTreshold;
            this.AgglomerationThreshold_Oldtimesteps = oldTs__AgglomerationTreshold;

            /*
            CutCellMetrics[] oldCcm;
            if (AgglomerateNewborn || AgglomerateDecased) {
                oldCcm = new CutCellMetrics[oldTs__AgglomerationTreshold.Length];
                for (int iHistory = 0; iHistory < oldCcm.Length; iHistory++) {
                    oldCcm[iHistory] = lsTrk.GetXDGSpaceMetrics(Spc, CutCellsQuadOrder, -iHistory).CutCellMetrics;
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

            /*
            {
                var grdDat = this.Tracker.GridDat;
                int Jup = grdDat.Cells.NoOfLocalUpdatedCells;

                Basis b = new Basis(grdDat, 0);
                var CellVolumesViz = new List<DGField>();
                for(int n = 0; n < oldCcm.Length; n++) {
                    foreach(var pair in oldCcm[n].CutCellVolumes) {
                        string name = this.Tracker.GetSpeciesName(pair.Key);
                        MultidimensionalArray vol = pair.Value;
                        var f = new SinglePhaseField(b, $"Volume#{name}AtIdx{n}");
                        CellVolumesViz.Add(f);

                        for(int j = 0; j < Jup; j++) {
                            f.SetMeanValue(j, vol[j]);
                        }



                    }
                 
                    var maskA = oldCcm[n].XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask("A");
                    var maskB = oldCcm[n].XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask("B");
                    var maskC = oldCcm[n].XDGSpaceMetrics.LevelSetRegions.GetSpeciesMask("C");

                    var fa = new SinglePhaseField(b, "mask-A"); fa.AccConstant(1.0, maskA);
                    var fb = new SinglePhaseField(b, "mask-B"); fb.AccConstant(1.0, maskB);
                    var fc = new SinglePhaseField(b, "mask-C"); fc.AccConstant(1.0, maskC);

                    CellVolumesViz.Add(fa);
                    CellVolumesViz.Add(fb);
                    CellVolumesViz.Add(fc);
                }

                Katastrophenplot(CellVolumesViz.ToArray());

                throw new Exception();
            }
            */

            // perform agglomeration
            foreach (var spc in this.SpeciesList) {
                /*
                IEnumerable<Tuple<int, int>> ai = FindAgglomeration(
                    this.Tracker,
                    spc,
                    AgglomerationThreshold,
                    this.NonAgglomeratedMetrics.CutCellVolumes[spc],
                    this.NonAgglomeratedMetrics.CutEdgeAreas[spc],
                    AgglomerateNewborn, AgglomerateDecased,
                    ExceptionOnFailedAgglomeration,
                    oldCcm != null ? oldCcm.Select(a => a.CutCellVolumes[spc]).ToArray() : null,
                    oldTs__AgglomerationTreshold,
                    NewbornAndDecasedThreshold);
                */
                                

                var aggAlg = new AgglomerationAlgorithm(this.Tracker, spc, CutCellsQuadOrder,
                    AgglomerationThreshold, oldTs__AgglomerationTreshold, NewbornAndDecasedThreshold,
                    AgglomerateNewborn, AgglomerateDecased,
                    ExceptionOnFailedAgglomeration
                    );

                var m_agglomeration = new CellAgglomerator(this.Tracker.GridDat, aggAlg.AgglomerationPairs);

                int myRank = lsTrk.GridDat.MpiRank;
                foreach(var p in m_agglomeration.AggInfo.AgglomerationPairs) {
                    Console.Error.WriteLine($"Rnk {myRank}, Spc {lsTrk.GetSpeciesName(spc)}: loc/gid {p.jCellSource}/{lsTrk.GridDat.iLogicalCells.GetGlobalID(p.jCellSource)} -> {p.jCellTarget}/{lsTrk.GridDat.iLogicalCells.GetGlobalID(p.jCellTarget)} (lv {p.AgglomerationLevel})");
                }


                this.DictAgglomeration.Add(spc, m_agglomeration);
            }

            this.TotalNumberOfAgglomerations = this.DictAgglomeration.Values.Sum(agg => agg.TotalNumberOfAgglomerations);

            // compute metrics of AGGLOMERATED cut cells
            this.LengthScaleAgg();
        }


        IDictionary<SpeciesId, CellAgglomerator> DictAgglomeration = new Dictionary<SpeciesId, CellAgglomerator>();


        /// <summary>
        /// returns the <see cref="CellAgglomerator"/> for species <paramref name="spc"/>;
        /// </summary>
        public CellAgglomerator GetAgglomerator(SpeciesId spc) {
            return DictAgglomeration[spc];
        }


        private Dictionary<SpeciesId, XSpatialOperatorMk2.SpeciesFrameMatrix<IMutableMatrixEx>> GetFrameMatrices<M>(M Matrix, UnsetteledCoordinateMapping RowMap, UnsetteledCoordinateMapping ColMap)
            where M : IMutableMatrixEx {
            var ret = new Dictionary<SpeciesId, XSpatialOperatorMk2.SpeciesFrameMatrix<IMutableMatrixEx>>();
            foreach (var kv in DictAgglomeration) {
                var Species = kv.Key;
                var mtx_spc = new XSpatialOperatorMk2.SpeciesFrameMatrix<IMutableMatrixEx>(Matrix, this.Tracker.Regions, Species, RowMap, ColMap);
                ret.Add(Species, mtx_spc);
            }

            return ret;
        }


        class MiniMapping {

            /// <summary>
            /// DOFs per cell, per variable, per species.
            /// index: variable.
            /// </summary>
            int[] NS;

            bool[] VarIsXdg;

            UnsetteledCoordinateMapping m_Map;
            SpeciesId m_spId;

            LevelSetTracker.LevelSetRegions m_LsRegion;

            public int MaxDeg = -1;

            public int NoOfVars;

            public MiniMapping(UnsetteledCoordinateMapping Map, SpeciesId spId, LevelSetTracker.LevelSetRegions r) {
                m_Map = Map;
                m_spId = spId;
                m_LsRegion = r;

                Basis[] BS = Map.BasisS.ToArray();
                NS = new int[BS.Length];
                VarIsXdg = new bool[BS.Length];
                NoOfVars = BS.Length;

                for (int iVar = 0; iVar < BS.Length; iVar++) {
                    XDGBasis xBasis = BS[iVar] as XDGBasis;
                    if (xBasis != null) {
                        NS[iVar] = xBasis.NonX_Basis.Length;
                        //m_LsTrk = xBasis.Tracker;
                        VarIsXdg[iVar] = true;
                    } else {
                        NS[iVar] = BS[iVar].Length;
                        VarIsXdg[iVar] = false;
                    }

                    MaxDeg = Math.Max(MaxDeg, BS[iVar].Degree);
                }
            }

            public long i0Func(int jCell, int iVar) {
                if (VarIsXdg[iVar]) {
                    int iSpc = m_LsRegion.GetSpeciesIndex(this.m_spId, jCell);
                    return m_Map.GlobalUniqueCoordinateIndex(iVar, jCell, iSpc * NS[iVar]);
                } else {
                    return m_Map.GlobalUniqueCoordinateIndex(iVar, jCell, 0);
                }
            }

            public int NFunc(int jCell, int iVar) {
                return NS[iVar];
            }

        }

        /// <summary>
        /// applies the agglomeration on a general matrix
        /// </summary>
        /// <param name="Matrix">the matrix that should be manipulated.</param>
        /// <param name="Rhs">the right-hand-side that should be manipulated</param>
        /// <param name="ColMap"></param>
        /// <param name="ColMapAggSw">Turns column agglomeration on/off fore each variable individually; default == null is on. </param>
        /// <param name="RowMap"></param>
        /// <param name="RowMapAggSw">The same shit as for <paramref name="ColMapAggSw"/>, just for rows.</param>
        public void ManipulateMatrixAndRHS<M, T>(M Matrix, T Rhs, UnsetteledCoordinateMapping RowMap, UnsetteledCoordinateMapping ColMap, bool[] RowMapAggSw = null, bool[] ColMapAggSw = null)
            where M : IMutableMatrixEx //
            where T : IList<double> //
        {
            using (var Ft = new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                //var mtxS = GetFrameMatrices(Matrix, RowMap, ColMap);

                Console.WriteLine("ManipulateMatrixAndRHS .................. ");

                if (Matrix == null && Rhs == null)
                    // nothing to do
                    return;

                if (TotalNumberOfAgglomerations <= 0)
                    // nothing to do
                    return;

                if (RowMapAggSw != null)
                    throw new NotImplementedException();

                // generate agglomeration sparse matrices
                // ======================================

                int RequireRight;
                if (Matrix == null) {
                    // we don't need multiplication-from-the-right at all
                    RequireRight = 0;
                } else {
                    if (RowMap.EqualsUnsetteled(ColMap) && ArrayTools.ListEquals(ColMapAggSw, RowMapAggSw)) {
                        // we can use the same matrix for right and left multiplication
                        RequireRight = 1;

                    } else {
                        // separate matrix for the multiplication-from-the-right is required
                        RequireRight = 2;
                    }
                }

                BlockMsrMatrix LeftMul = null, RightMul = null;
                {

                    foreach (var kv in DictAgglomeration) {
                        var Species = kv.Key;
                        var m_Agglomerator = kv.Value;

                        if (m_Agglomerator != null) {

                            CellMask spcMask = this.Tracker.Regions.GetSpeciesMask(Species);

                            MiniMapping rowMini = new MiniMapping(RowMap, Species, this.Tracker.Regions);
                            BlockMsrMatrix LeftMul_Species = m_Agglomerator.GetRowManipulationMatrix(RowMap, rowMini.MaxDeg, rowMini.NoOfVars, rowMini.i0Func, rowMini.NFunc, false, spcMask, this.Tracker.GetSpeciesName(Species));
                            if (LeftMul == null) {
                                LeftMul = LeftMul_Species;
                            } else {
                                LeftMul.Acc(1.0, LeftMul_Species);
                            }


                            if (!object.ReferenceEquals(LeftMul, RightMul) && RightMul != null) {
                                MiniMapping colMini = new MiniMapping(ColMap, Species, this.Tracker.Regions);
                                BlockMsrMatrix RightMul_Species = m_Agglomerator.GetRowManipulationMatrix(ColMap, colMini.MaxDeg, colMini.NoOfVars, colMini.i0Func, colMini.NFunc, false, spcMask, this.Tracker.GetSpeciesName(Species));

                                if (RightMul == null) {
                                    RightMul = RightMul_Species;
                                } else {
                                    RightMul.Acc(1.0, RightMul_Species);
                                }

                            } else if (RequireRight == 1) {
                                RightMul = LeftMul;
                            } else {
                                RightMul = null;
                            }
                        }
                    }
                }

                // apply the agglomeration to the matrix
                // =====================================

                if (Matrix != null) {
                    BlockMsrMatrix RightMulTr = RightMul.Transpose();

                    BlockMsrMatrix _Matrix;
                    if (Matrix is BlockMsrMatrix) {
                        _Matrix = (BlockMsrMatrix)((object)Matrix);
                    } else {
                        _Matrix = Matrix.ToBlockMsrMatrix(RowMap, ColMap);
                    }

                    var AggMatrix = BlockMsrMatrix.Multiply(LeftMul, BlockMsrMatrix.Multiply(_Matrix, RightMulTr));

                    if (object.ReferenceEquals(_Matrix, Matrix)) {
                        _Matrix.Clear();
                        _Matrix.Acc(1.0, AggMatrix);
                    } else {
                        Matrix.Acc(-1.0, _Matrix); //   das ist so
                        Matrix.Acc(1.0, AggMatrix); //  meagaschlecht !!!!!!
                    }
                }

                // apply the agglomeration to the Rhs
                // ==================================

                if (Rhs != null) {

                    double[] tmp = Rhs.ToArray();
                    if (object.ReferenceEquals(tmp, Rhs))
                        throw new ApplicationException("Flache kopie sollte eigentlich ausgeschlossen sein!?");

                    LeftMul.SpMV(1.0, tmp, 0.0, Rhs);
                }

                Console.WriteLine("........................ ManipulateMatrixAndRHS ");
            }
        }

        Dictionary<SpeciesId, XSpatialOperatorMk2.SpeciesFrameVector<T>> GetFrameVectors<T>(T vec, UnsetteledCoordinateMapping Map)
            where T : IList<double> {
            var ret = new Dictionary<SpeciesId, XSpatialOperatorMk2.SpeciesFrameVector<T>>();
            foreach (var kv in DictAgglomeration) {
                var Species = kv.Key;
                var vec_spc = new XSpatialOperatorMk2.SpeciesFrameVector<T>(this.Tracker.Regions, Species, vec, Map);
                ret.Add(Species, vec_spc);
            }
            return ret;
        }


        /// <summary>
        /// for all DG fields in <paramref name="Map"/>, this clears all entries which correspond to agglomerated cells.
        /// </summary>
        public void ClearAgglomerated(CoordinateMapping Map) {
            var vec = new CoordinateVector(false, Map.Fields.ToArray());
            ClearAgglomerated(vec, Map);
        }

        /// <summary>
        /// In a vector <paramref name="vec"/>, this clears all entries which correspond to agglomerated cells.
        /// </summary>
        public void ClearAgglomerated<T>(T vec, UnsetteledCoordinateMapping Map)
            where T : IList<double> {
            var vecS = GetFrameVectors(vec, Map);

            foreach (var kv in DictAgglomeration) {
                var Species = kv.Key;
                var m_Agglomerator = kv.Value;


                if (m_Agglomerator != null) {
                    var vec_spc = vecS[Species];
                    m_Agglomerator.ClearAgglomerated(vec_spc, vec_spc.Mapping);
                }

            }
        }

        /// <summary>
        /// computes the l2-Norm of the DG coordinates
        /// in the agglomerated cells
        /// </summary>
        /// <remarks>
        /// The l2-Norm of the DG coordinates
        /// is NOT equal to the L2-Norm of the DG-field in the XDG-case,
        /// since the basis is no longer orthonormal and therefore Parcival's equation does not hold.
        /// </remarks>
        public double NormInAgglomerated<T>(T vec, UnsetteledCoordinateMapping Map)
            where T : IList<double> {

            var vecS = GetFrameVectors(vec, Map);

            double l2Norm = 0.0;
            foreach (var kv in DictAgglomeration) {
                var Species = kv.Key;
                var m_Agglomerator = kv.Value;


                if (m_Agglomerator != null) {
                    var vec_spc = vecS[Species];
                    l2Norm += m_Agglomerator.NormInAgglomerated(vec_spc, vec_spc.Mapping).Pow2();
                }
            }

            return l2Norm.Sqrt();
        }

        /// <summary>
        /// In a vector <paramref name="vec"/>, this method performs a
        /// polynomial extrapolation from agglomeration target cells to agglomeration source cells.
        /// </summary>
        public void Extrapolate(CoordinateMapping Map) {
            //var vecS = GetFrameVectors(vec, Map);

            foreach (var kv in DictAgglomeration) {
                var Species = kv.Key;
                var m_Agglomerator = kv.Value;

                var DgFields = Map.Fields.ToArray();

                DGField[] SubFields = new DGField[DgFields.Count()];
                for (int iFld = 0; iFld < SubFields.Length; iFld++) {
                    DGField f = DgFields[iFld];

                    if (f is ConventionalDGField)
                        SubFields[iFld] = (ConventionalDGField)(f);
                    else if (f is XDGField)
                        SubFields[iFld] = ((XDGField)f).GetSpeciesShadowField(Species);
                    else
                        throw new NotImplementedException();
                }

                if (m_Agglomerator != null) {
                    m_Agglomerator.Extrapolate(new CoordinateMapping(SubFields));
                }

            }
        }


        /// <summary>
        /// Diagnostic output: for each species, 'stays' which visualize the agglomeration operation, 
        /// are plotted.
        /// </summary>
        /// <param name="basename">base name for the output files</param>
        public void PlotAgglomerationPairs(string basename) {

            var xqs = this.XDGSpaceMetrics.XQuadSchemeHelper;

            foreach (var Species in this.Tracker.SpeciesIdS) {
                if (!DictAgglomeration.ContainsKey(Species)) {
                    continue;
                }

                var gdat = this.XDGSpaceMetrics.GridDat;
                int J = gdat.Cells.NoOfLocalUpdatedCells;
                int JE = gdat.Cells.Count;
                int D = gdat.SpatialDimension;

                var CompScheme = xqs.GetVolumeQuadScheme(Species).Compile(gdat, this.CutCellQuadratureOrder);

                MultidimensionalArray CenterOfGravity = MultidimensionalArray.Create(JE, D);

                BoSSS.Foundation.Quadrature.CellQuadrature.GetQuadrature(new int[] { D + 1 }, gdat,
                    CompScheme,
                    delegate (int i0, int Length, QuadRule rule, MultidimensionalArray EvalResult) { // Del_Evaluate
                        var globNodes = gdat.GlobalNodes.GetValue_Cell(rule.Nodes, i0, Length);
                        EvalResult.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Length - 1, rule.NoOfNodes - 1, D - 1 }).Set(globNodes);
                        EvalResult.ExtractSubArrayShallow(new int[] { 0, 0, D }, new int[] { Length - 1, rule.NoOfNodes - 1, D - 1 }).SetAll(1.0);
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++) {
                            for (int d = 0; d < D; d++) {
                                CenterOfGravity[i + i0, d] = ResultsOfIntegration[i, d] / ResultsOfIntegration[i, D];
                            }
                        }
                    }, cs: CoordinateSystem.Physical).Execute();

                CenterOfGravity.Storage.MPIExchange(this.Tracker.GridDat);


                var SpeciesName = this.Tracker.GetSpeciesName(Species);
                this.DictAgglomeration[Species].PlotAgglomerationPairs(SpeciesName + "-" + basename + ".csv", CenterOfGravity);
            }
        }

        /// <summary>
        /// The volume over cut cell surface ratio, i.e. \f$ \frac{ | K^X |}{ | \partial K^X | } \f$, for each **agglomerated** cut-cell $K^X$.
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CellLengthScales {
            private set;
            get;
        }

        /// <summary>
        /// The volume fraction of agglomerated cut cells.
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CellVolumeFrac {
            private set;
            get;
        }

        /// <summary>
        /// The volume of agglomerated cut cells.
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CutCellVolumes {
            private set;
            get;
        }

        /// <summary>
        /// The surface of agglomerated cut cells.
        /// </summary>
        public Dictionary<SpeciesId, MultidimensionalArray> CellSurface {
            private set;
            get;
        }

        /// <summary>
        /// Initializes <see cref="CellLengthScales"/>.
        /// </summary>
        void LengthScaleAgg() {
            using (new FuncTrace()) {
                SpeciesId[] species = this.SpeciesList.ToArray();

                int J = this.Tracker.GridDat.Cells.NoOfLocalUpdatedCells;
                int JE = this.Tracker.GridDat.Cells.Count;
                int[][] C2E = this.Tracker.GridDat.Cells.Cells2Edges;

                //TestingIO Checker = CheckFile != null ? new TestingIO(this.Tracker.GridDat, CheckFile, 1) : null;


                var CellLengthScalesMda = MultidimensionalArray.Create(JE, species.Length, 2); // 1st index: cell, 2nd index: species, 3rd index: [surface, volume]
                var CellVolumeFracMda = MultidimensionalArray.Create(JE, species.Length); // 1st index: cell, 2nd index: species

                for (int iSpc = 0; iSpc < species.Length; iSpc++) {
                    SpeciesId spc = species[iSpc];
                    var agginfo = this.GetAgglomerator(spc).AggInfo;
                    BitArray aggEdgesBitMask = agginfo.AgglomerationEdges.GetBitMask();

                    MultidimensionalArray CellSurface = CellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc, 0);
                    MultidimensionalArray CellVolume = CellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc, 1);

                    MultidimensionalArray CellVolume2 = CellVolumeFracMda.ExtractSubArrayShallow(-1, iSpc);

                    CellSurface.Set(this.NonAgglomeratedMetrics.InterfaceArea[spc]);
                    CellVolume.Set(this.NonAgglomeratedMetrics.CutCellVolumes[spc]);
                    CellVolume2.Set(this.NonAgglomeratedMetrics.CutCellVolumes[spc]);



                    MultidimensionalArray EdgeArea = this.NonAgglomeratedMetrics.CutEdgeAreas[spc];

                    // accumulate cell surface
                    for (int j = 0; j < J; j++) {
                        int[] edges = C2E[j];
                        int NE = edges.Length;

                        for (int ne = 0; ne < NE; ne++) {
                            int iEdg = Math.Abs(edges[ne]) - 1;

                            Debug.Assert(this.Tracker.GridDat.Edges.CellIndices[iEdg, 0] == j || this.Tracker.GridDat.Edges.CellIndices[iEdg, 1] == j);

                            if (!aggEdgesBitMask[iEdg]) { // exclude edges in agglomeration pairs
                                Debug.Assert(!(double.IsNaN(CellSurface[j]) || double.IsInfinity(CellSurface[j])));
                                Debug.Assert(!(double.IsNaN(EdgeArea[iEdg]) || double.IsInfinity(EdgeArea[iEdg])));
                                CellSurface[j] += EdgeArea[iEdg];
                                Debug.Assert(!(double.IsNaN(CellSurface[j]) || double.IsInfinity(CellSurface[j])));
                            }
                        }
                    }

                    //if(Checker != null) {
                    //    Checker.AddVector("CellSurface" + this.Tracker.GetSpeciesName(spc), CellSurface.To1DArray().GetSubVector(0, J));
                    //    Checker.AddVector("CellVolume" + this.Tracker.GetSpeciesName(spc), CellVolume.To1DArray().GetSubVector(0, J));
                    //    Checker.AddVector("CellVolume2" + this.Tracker.GetSpeciesName(spc), CellVolume2.To1DArray().GetSubVector(0, J));
                    //}
                }

                //if(Checker != null) {
                //    Checker.DoIOnow();
                //    foreach(string cn in Checker.ColumnNamesWithoutReserved) {
                //        double d = Checker.RelError(cn);
                //        Console.WriteLine($"   ------------ rel error of {cn} before comm: " + d);
                //    }

                //    Checker.CurrentData.Clear();
                //}

                // MPI exchange
                // Needed, such that all ExternalCells (i.e. Ghost cells) have the correct CellSurface
                CellLengthScalesMda.Storage.MPIExchange(this.Tracker.GridDat);
                CellVolumeFracMda.Storage.MPIExchange(this.Tracker.GridDat);

                //if(Checker != null) {
                //    for(int iSpc = 0; iSpc < species.Length; iSpc++) {
                //        SpeciesId spc = species[iSpc];

                //        MultidimensionalArray CellSurface = CellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc, 0);
                //        MultidimensionalArray CellVolume = CellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc, 1);

                //        MultidimensionalArray CellVolume2 = CellVolumeFracMda.ExtractSubArrayShallow(-1, iSpc);


                //        Checker.AddVector("CellSurface" + this.Tracker.GetSpeciesName(spc), CellSurface.To1DArray().GetSubVector(0, J));
                //        Checker.AddVector("CellVolume" + this.Tracker.GetSpeciesName(spc), CellVolume.To1DArray().GetSubVector(0, J));
                //        Checker.AddVector("CellVolume2" + this.Tracker.GetSpeciesName(spc), CellVolume2.To1DArray().GetSubVector(0, J));
                //    }

                //    foreach(string cn in Checker.ColumnNamesWithoutReserved) {
                //        double d = Checker.RelError(cn);
                //        Console.WriteLine($"   ------------ rel error of {cn} AFTER comm: " + d);
                //    }
                //}

                var AggCellLengthScalesMda = MultidimensionalArray.Create(JE, species.Length); // 1st index: cell, 2nd index: species
                for (int iSpc = 0; iSpc < species.Length; iSpc++) {
                    SpeciesId spc = species[iSpc];
                    var agginfo = this.GetAgglomerator(spc).AggInfo;

                    MultidimensionalArray CellSurface = CellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc, 0);
                    MultidimensionalArray CellVolume = CellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc, 1);

                    MultidimensionalArray CellVolume2 = CellVolumeFracMda.ExtractSubArrayShallow(-1, iSpc);

                    // sum agglomeration sources to targets
                    foreach (var agg_pair in agginfo.AgglomerationPairs) {
                        CellSurface[agg_pair.jCellTarget] += CellSurface[agg_pair.jCellSource];
                        CellVolume[agg_pair.jCellTarget] += CellVolume[agg_pair.jCellSource];
                    }

                    MultidimensionalArray LengthScales = AggCellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc);

                    // Loop includes external cells
                    for (int j = 0; j < JE; j++) {
                        LengthScales[j] = CellVolume[j] / CellSurface[j];
                        CellVolume2[j] = CellVolume2[j] / this.Tracker.GridDat.Cells.GetCellVolume(j);
                    }

                    // set values in agglomeration sources to be equal to agglomeration targets
                    foreach (var agg_pair in agginfo.AgglomerationPairs) {
                        LengthScales[agg_pair.jCellSource] = LengthScales[agg_pair.jCellTarget];
                        CellVolume2[agg_pair.jCellSource] = CellVolume2[agg_pair.jCellTarget];
                    }

                    if (this.AgglomerationThreshold <= 0.0) {
                        // special treatment for no agglomeration -- which is anyway not recommended at all
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        CellMask spcDom = this.Tracker.Regions.GetSpeciesMask(spc);

                        foreach (int j in spcDom.ItemEnum) {

                            double unCut = this.Tracker.GridDat.Cells.GetCellVolume(j);
                            double Fraction = CellVolume[j] / unCut;

                            if (Fraction <= 1.0e-10)
                                LengthScales[j] = 1e10;
                        }
                    }


                }

                // MPI exchange -> Is it really needed now???
                // Yes! we need length scales for external/ghost cells, in order to compute fluxes at the boundaries
                AggCellLengthScalesMda.Storage.MPIExchange(this.Tracker.GridDat);

                // store
                this.CellLengthScales = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.CellVolumeFrac = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.CellSurface = new Dictionary<SpeciesId, MultidimensionalArray>();
                this.CutCellVolumes = new Dictionary<SpeciesId, MultidimensionalArray>();
                var LsChecker = new TestingIO(this.Tracker.GridDat, $"LevelSets-LenScale.abc", 1);
                for (int iSpc = 0; iSpc < species.Length; iSpc++) {
                    SpeciesId spc = species[iSpc];
                    this.CellLengthScales.Add(spc, AggCellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc).CloneAs());
                    this.CellVolumeFrac.Add(spc, CellVolumeFracMda.ExtractSubArrayShallow(-1, iSpc).CloneAs());
                    this.CellSurface.Add(spc, CellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc, 0).CloneAs());
                    this.CutCellVolumes.Add(spc, CellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc, 1).CloneAs());


                    //for(int j = 0; j < J; j++) {
                    //    Console.Error.WriteLine($"Rnk {this.Tracker.GridDat.MpiRank}, Spc {this.Tracker.GetSpeciesName(spc)} gid {this.Tracker.GridDat.iLogicalCells.GetGlobalID(j)}: {AggCellLengthScalesMda[j, iSpc]} = {CellLengthScalesMda[j, iSpc, 0]} {CellLengthScalesMda[j, iSpc, 1]} ");
                    //}

                    LsChecker.AddVector("LenScale-" + this.Tracker.GetSpeciesName(spc), this.CellLengthScales[spc].To1DArray().Take(J).Select(a => a.IsNaN() ? -99.1 : a));
                    LsChecker.AddVector("Vol-" + this.Tracker.GetSpeciesName(spc), this.CutCellVolumes[spc].To1DArray().Take(J).Select(a => a.IsNaN() ? -99.2 : a));
                    LsChecker.AddVector("Surf-" + this.Tracker.GetSpeciesName(spc), this.CellSurface[spc].To1DArray().Take(J).Select(a => a.IsNaN() ? -99.3 : a));

                    this.CellLengthScales[spc].SetAll(0.1);
                    this.CellVolumeFrac[spc].SetAll(0.1);
                    this.CellSurface[spc].SetAll(0.1);
                    this.CutCellVolumes[spc].SetAll(0.1);
                }

                LsChecker.DoIOnow();
                var err = LsChecker.AllAbsErr();
                foreach(var kv in err) {
                    Console.WriteLine($"    Err {kv.Key} = {kv.Value}");
                }
            }
        }

        /*
        /// <summary>
        /// Sometimes, this provides indeed a correct agglomeration graph.
        /// If agglomeration fails -- which it does quite regularly -- unleash the <see cref="Katastrophenplot"/> to see the mess!
        /// </summary>
        /// <remarks>
        /// Cell agglomeration is used to handle two problems:
        /// first, for the treatment of very small cut cells and, for temporally evolving interfaces, to 
        /// ensure an equal topology of the (agglomerated) XDG cut-cell mesh for all involved temporal levels.
        /// - the issue of small cut cells is described in the paper:
        ///   _Extended discontinuous Galerkin methods for two-phase flows: the spatial discretization; Kummer; IJNMF 109 (2), 2017_. 
        /// - the agglomeration of _newborn_ and _decased_ cells is described in 
        ///   the paper: _Time integration for extended discontinuous Galerkin methods with moving domains; Kummer, Müller, Utz; IJNMF 113 (5), 2018_.
        /// </remarks>
        static public IEnumerable<Tuple<int, int>> FindAgglomeration(LevelSetTracker Tracker, SpeciesId spId, double AgglomerationThreshold,
            MultidimensionalArray CellVolumes, MultidimensionalArray edgeArea,
            bool AgglomerateNewborn, bool AgglomerateDeceased, bool ExceptionOnFailedAgglomeration,
            MultidimensionalArray[] oldCellVolumes, double[] oldTs__AgglomerationTreshold, double NewbornAndDecasedThreshold) //
        {

            var as_data = FindAgglomerationSources(Tracker, spId, AgglomerationThreshold,
                CellVolumes, 
                AgglomerateNewborn, AgglomerateDeceased,
                oldCellVolumes, oldTs__AgglomerationTreshold, NewbornAndDecasedThreshold);

            var ret = FindAgglomerationTargets(Tracker, spId, CellVolumes, edgeArea, oldCellVolumes, as_data.AgglomCellsList, as_data.AgglomCellsBitmask, as_data.AggCandidates, ExceptionOnFailedAgglomeration);

            return ret;
        }
        */

       
    }
}
