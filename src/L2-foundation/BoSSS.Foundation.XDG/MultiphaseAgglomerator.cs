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
        /// Agglomeration thresholds for   _previous_ timesteps, correlates with <paramref name="oldCcm"/>.
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
        /// <param name="NewbornAndDecasedThreshold">
        /// Volume fraction threshold at which a cut-cell counts as newborn, resp. deceased, see <paramref name="AgglomerateNewbornAndDeceased"/>;
        /// </param>
        /// <param name="CutCellsQuadOrder"></param>
        /// <param name="Spc"></param>
        /// <param name="lsTrk"></param>
        internal MultiphaseCellAgglomerator(
            LevelSetTracker lsTrk,
            SpeciesId[] Spc, int CutCellsQuadOrder,
            double __AgglomerationTreshold,
            bool AgglomerateNewborn = false, bool AgglomerateDecased = false, bool ExceptionOnFailedAgglomeration = true,
            double[] oldTs__AgglomerationTreshold = null,
            double NewbornAndDecasedThreshold = 1.0e-6) {
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


                var m_agglomeration = new CellAgglomerator(this.Tracker.GridDat, ai);
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
                            BlockMsrMatrix LeftMul_Species = m_Agglomerator.GetRowManipulationMatrix(RowMap, rowMini.MaxDeg, rowMini.NoOfVars, rowMini.i0Func, rowMini.NFunc, false, spcMask);
                            if (LeftMul == null) {
                                LeftMul = LeftMul_Species;
                            } else {
                                LeftMul.Acc(1.0, LeftMul_Species);
                            }


                            if (!object.ReferenceEquals(LeftMul, RightMul) && RightMul != null) {
                                MiniMapping colMini = new MiniMapping(ColMap, Species, this.Tracker.Regions);
                                BlockMsrMatrix RightMul_Species = m_Agglomerator.GetRowManipulationMatrix(ColMap, colMini.MaxDeg, colMini.NoOfVars, colMini.i0Func, colMini.NFunc, false, spcMask);

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
        /// The volume over cut cell surface ratio, i.e. \f$ \frac{ | K^X |}{ | \partial K^X | } \f$, for each agglomerated cut-cell $K^X$.
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
                for (int iSpc = 0; iSpc < species.Length; iSpc++) {
                    SpeciesId spc = species[iSpc];
                    this.CellLengthScales.Add(spc, AggCellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc).CloneAs());
                    this.CellVolumeFrac.Add(spc, CellVolumeFracMda.ExtractSubArrayShallow(-1, iSpc).CloneAs());
                    this.CellSurface.Add(spc, CellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc, 0).CloneAs());
                    this.CutCellVolumes.Add(spc, CellLengthScalesMda.ExtractSubArrayShallow(-1, iSpc, 1).CloneAs());
                }
            }
        }

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

        static (List<int> AgglomCellsList, BitArray AgglomCellsBitmask, BitArray AggCandidates) FindAgglomerationSources(LevelSetTracker Tracker, SpeciesId spId, double AgglomerationThreshold,
            MultidimensionalArray CellVolumes,
            bool AgglomerateNewborn, bool AgglomerateDeceased,
            MultidimensionalArray[] oldCellVolumes, double[] oldTs__AgglomerationTreshold, double NewbornAndDecasedThreshold) //
        {

            using (var tracer = new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                tracer.InfoToConsole = true;
                tracer.Info("agglom newborn? " + AgglomerateNewborn);
                tracer.Info("agglom decased? " + AgglomerateDeceased);
                tracer.Info("AgglomerationThreshold = " + AgglomerationThreshold + ",  NewbornAndDecasedThreshold = " + NewbornAndDecasedThreshold);
                if (oldTs__AgglomerationTreshold != null)
                    tracer.Info("oldTs__AgglomerationTreshold = " + oldTs__AgglomerationTreshold.ToConcatString("", "| ", ";"));

                // init 
                // =======

                GridData grdDat = Tracker.GridDat;
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
                            Console.WriteLine("Must agglom cell " + jCell + "#" + Tracker.GetSpeciesName(spId) + " volume frac is " + frac);
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

                                    Console.WriteLine("Must agglom cell " + jCell + "#" + Tracker.GetSpeciesName(spId) + " volume frac is " + frac + "on time level " + iTimeLev);

                                }
                            }
                        }
                    }
                }

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
                                    AggCandidates[jNewbornCell] = false;
                                    if (!AgglomCellsBitmask[jNewbornCell]) {
                                        AgglomCellsList.Add(jNewbornCell);
                                        AgglomCellsBitmask[jNewbornCell] = true;

                                        Console.WriteLine("Must agglom NEWBORN cell " + AgglomCellsBitmask + "#" + Tracker.GetSpeciesName(spId));

                                    }
                                }
                            }
                        }
                    }
                }

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

                                    int jNewbornCell = j;
                                    AggCandidates[jNewbornCell] = false;
                                    if (!AgglomCellsBitmask[jNewbornCell]) {
                                        AgglomCellsList.Add(jNewbornCell);
                                        AgglomCellsBitmask[jNewbornCell] = true;

                                        Console.WriteLine("Must agglom DEAD cell " + AgglomCellsBitmask + "#" + Tracker.GetSpeciesName(spId));

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

                return (AgglomCellsList, AgglomCellsBitmask, AggCandidates);

            }
        }

        static IEnumerable<Tuple<int, int>> FindAgglomerationTargets(LevelSetTracker Tracker, SpeciesId spId,
            MultidimensionalArray CellVolumes, MultidimensionalArray edgeArea, MultidimensionalArray[] oldCellVolumes,
            List<int> AgglomCellsList, BitArray AgglomCellsBitmask, BitArray AggCandidates, bool ExceptionOnFailedAgglomeration
            ) {
            using (new FuncTrace()) {

                GridData grdDat = Tracker.GridDat;
                var Cell2Edge = grdDat.Cells.Cells2Edges;
                int[,] Edge2Cell = grdDat.Edges.CellIndices;
                int NoOfEdges = grdDat.Edges.Count;
                byte[] EdgeTags = grdDat.Edges.EdgeTags;
                int myMpiRank = Tracker.GridDat.MpiRank;
                int Jup = grdDat.Cells.NoOfLocalUpdatedCells;
                int Jtot = grdDat.Cells.Count;
                Partitioning CellPart = Tracker.GridDat.CellPartitioning;
                var GidxExt = Tracker.GridDat.Parallel.GlobalIndicesExternalCells;


                var AgglomerationPairs = new List<CellAgglomerator.AgglomerationPair>();

                if (edgeArea.GetLength(0) != NoOfEdges)
                    throw new ArgumentException();

                double EmptyEdgeTreshold = 1.0e-10; // edges with a measure below or equal to this threshold are
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

                    bool print = false;
                    if (jCell == 29 || jCell == 30) {
                        print = true;
                        Console.WriteLine("Looking for agglom for cell " + jCell + "#" + Tracker.GetSpeciesName(spId) + " volume is " + CellVolumes[jCell]);
                    }


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
                        if (print) {
                            Console.WriteLine("  testing with cell " + jCellNeigh);
                            Console.WriteLine("    connecting edge area: " + EdgeArea_iEdge);
                        }
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

                        if (print)
                            Console.WriteLine("    neighbour fraction: " + frac_neigh);

                        // max?
                        if (frac_neigh > frac_neigh_max) {
                            frac_neigh_max = frac_neigh;
                            e_max = e;
                            jCellNeigh_max = jCellNeigh;
                            jEdge_max = iEdge;
                        }
                    }

                    if (print)
                        Console.WriteLine("  selected: " + jCellNeigh_max);

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

                if (failCells.Count.MPISum() > 0) {
                    // ++++++++++++++++++++++++++
                    // Error handling / reporting
                    // ++++++++++++++++++++++++++

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
                    foreach (int j in failCells) {
                        FailedViz.SetMeanValue(j, 1);
                    }

                    DGField[] LevelSets = Tracker.LevelSets.Select(s => (DGField)s).ToArray();

                    if (Katastrophenplot != null)
                        Katastrophenplot(CellVolumesViz.Cat(AgglomCellsViz, FailedViz, LevelSets));



                    string message = ("Agglomeration failed - no candidate for agglomeration found");
                    if (ExceptionOnFailedAgglomeration)
                        throw new Exception(message);
                    else
                        Console.WriteLine(message);



                }



                // store & return
                // ================
                return AgglomerationPairs.Select(pair => new Tuple<int, int>(pair.jCellSource, pair.jCellTarget)).ToArray();

            }
        }
    

        /// <summary>
        /// Temporary feature; will be removed in future;
        /// Plotting if agglomeration fails.
        /// </summary>
        public static Action<DGField[]> Katastrophenplot;
    }
}
