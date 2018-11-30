using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using System.Diagnostics;
using ilPSP.Tracing;
using ilPSP.Utils;

namespace BoSSS.Foundation.Grid.Aggregation {
    public partial class AggregationGridData {
        static public Action<Classic.GridCommons, int[]> PlotScheisse;

        public BasisData ChefBasis {
            get {
                return m_ChefBasis;
            }
        }

        _BasisData m_ChefBasis;

        class _BasisData : BasisData {

            AggregationGridData m_owner;

            internal _BasisData(AggregationGridData o) : base(o) {
                m_owner = o;
            }

            int PolyDim(int Degree) {
                var Krefs = m_owner.m_GeomCellData.RefElements;
                var Polys = Krefs[0].GetOrthonormalPolynomials(Degree);
                return Polys.Count;
            }


            MultidimensionalArray m_OrthonormalizationTrafo;
            int m_OrthonormalizationTrafo_Degree = -1;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="j0"></param>
            /// <param name="Len"></param>
            /// <param name="Degree"></param>
            /// <returns></returns>
            protected override MultidimensionalArray Compute_OrthonormalizationTrafo(int j0, int Len, int Degree) {
                this.Init(Degree);

                int Np = PolyDim(Degree);

                if (m_OrthonormalizationTrafo == null || m_OrthonormalizationTrafo_Degree < Degree) {
                    // ++++++++++++++++++++++++++
                    // require re-evaluation
                    // ++++++++++++++++++++++++++

                    int Jgeom = m_owner.m_GeomCellData.Count;
                    m_OrthonormalizationTrafo = MultidimensionalArray.Create(Jgeom, Np, Np);

                    int Jlog = m_owner.m_LogicalCellData.Count;

                    // logical to geometrical transformation
                    int[][] jL2jG = m_owner.iLogicalCells.AggregateCellToParts;

                    for (int jlog = 0; jlog < Jlog; jlog++) {
                        var Trafo = this.CA(jlog, Np);

                        int[] jGS = jL2jG[jlog];
                        Debug.Assert(jGS.Length == Trafo.GetLength(0));

                        for (int i = 0; i < jGS.Length; i++) {
                            int jG = jGS[i];

                            m_OrthonormalizationTrafo.ExtractSubArrayShallow(jG, -1, -1)
                                .Acc(1.0, Trafo.ExtractSubArrayShallow(i, -1, -1));
                        }
                    }

                }

                return m_OrthonormalizationTrafo.ExtractSubArrayShallow(new int[] { j0, 0, 0 }, new int[] { j0 + Len - 1, Np - 1, Np - 1 });
            }

            MultidimensionalArray InjectorsBase;

            /// <summary>
            /// injectors to upper level
            /// - outer array index: aggregation cell index
            /// For each multidimensional array:
            /// - 1st index: enumeration of parts of upper level
            /// - 2nd index: row index
            /// - 3rd index: column index
            /// </summary>
            public MultidimensionalArray[] Injectors;

            /// <summary>
            /// Defines which entries of <see cref=InjectorsBase""/> contain valid values.
            /// - index: correlates with 1st index of <see cref="InjectorsBase"/>
            /// </summary>
            bool[] InjectorsBaseReady;

            int MaxSupportedDegree = -1;


            

            public void Init(int ReqDegree) {
                if (ReqDegree < 0)
                    throw new ArgumentOutOfRangeException();
                if (ReqDegree <= MaxSupportedDegree)
                    return;

                var aGdat = m_owner.AncestorGrid;
                Basis maxDGbasis = new Basis(aGdat, ReqDegree);
                int Np = maxDGbasis.Length;
                int Jbase = aGdat.Cells.NoOfLocalUpdatedCells;

                int Jagg = m_owner.iLogicalCells.NoOfLocalUpdatedCells;
                int[][] Ag2Pt = m_owner.iLogicalCells.AggregateCellToParts;
                int[][] C2F = m_owner.jCellCoarse2jCellFine;

                if (m_owner.ParentGrid is AggregationGridData) {
                    var agParrent = m_owner.ParentGrid as AggregationGridData;
                    int[][] Ag2Pt_Fine = agParrent.iLogicalCells.AggregateCellToParts;

                    agParrent.m_ChefBasis.Init(MaxSupportedDegree);

                    if (ReqDegree == agParrent.m_ChefBasis.MaxSupportedDegree) {
                        this.InjectorsBase = agParrent.m_ChefBasis.InjectorsBase.CloneAs();
                    } else {
                        this.InjectorsBase = agParrent.m_ChefBasis.InjectorsBase.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Jbase - 1, Np - 1, Np - 1 }).CloneAs();
                    }
                    this.InjectorsBaseReady = agParrent.m_ChefBasis.InjectorsBaseReady.CloneAs();

                    this.Injectors = BuildInjector_Lv2andup(maxDGbasis, Np, this.InjectorsBase, this.InjectorsBaseReady, Jagg, Ag2Pt_Fine, C2F);

                } else if (m_owner.ParentGrid is Classic.GridData) {
                    // this is level 0

                    var agParrent = m_owner.ParentGrid as Classic.GridData;
                    this.InjectorsBase = MultidimensionalArray.Create(Jbase, Np, Np);
                    this.InjectorsBaseReady = new bool[Jbase];
                    ArrayTools.SetAll(this.InjectorsBaseReady, true);

                    for (int j = 0; j < Jbase; j++) {
                        this.InjectorsBase.ExtractSubArrayShallow(j, -1, -1).AccEye(1.0);
                    }

                    if (IsCloneOfAncestor()) {
                        //this.Injectors = BuildInjector_Lv1(maxDGbasis, Np, this.InjectorsBase, this.InjectorsBaseReady, Jagg, Ag2Pt, C2F);
                        // nop
                    } else {
                        this.Injectors = BuildInjector_Lv1(maxDGbasis, Np, this.InjectorsBase, this.InjectorsBaseReady, Jagg, Ag2Pt, C2F);
                        //this.Injectors = BuildInjector_Lv2andup(maxDGbasis, Np, this.InjectorsBase, this.InjectorsBaseReady, Jagg, Ag2Pt_Fine, C2F);
                    }
                } else {
                    throw new NotSupportedException("dont know what to do");
                }

                MaxSupportedDegree = ReqDegree;
            }

            bool IsCloneOfAncestor() {
                var aGdat = m_owner.AncestorGrid;

                int Jbase = aGdat.Cells.NoOfLocalUpdatedCells;

                int Jagg = m_owner.iLogicalCells.NoOfLocalUpdatedCells;
                int[][] Ag2Pt = m_owner.iLogicalCells.AggregateCellToParts;
                if (Jagg != Jbase)
                    return false;

                for (int j = 0; j < Jagg; j++) {
                    int jGeom;
                    if (Ag2Pt == null || Ag2Pt[j] == null) {
                        jGeom = j;
                    } else {
                        if (Ag2Pt[j].Length != 1)
                            return false;
                        jGeom = Ag2Pt[j][0];
                    }
                    if (jGeom != j)
                        return false;
                }

                return true;
            }




            /// <summary>
            /// computes the injector for multigrid level 1
            /// </summary>
            /// <seealso cref="BuildInjector_Lv2andup"/>
            private static MultidimensionalArray[] BuildInjector_Lv1(
                Basis maxDgBasis, int Np,
                MultidimensionalArray InjectorsBase, bool[] InjectorsBaseReady,
                int Jagg, int[][] Ag2Pt, int[][] C2F) {
                using (new FuncTrace()) {
                    MultidimensionalArray ortho = MultidimensionalArray.Create(Np, Np);

                    MultidimensionalArray[] Injectors_iLevel = new MultidimensionalArray[Jagg];

#if DEBUG
                    {
                        int Jbase = InjectorsBase.GetLength(0);
                        if (InjectorsBase.GetLength(1) != InjectorsBase.GetLength(2))
                            throw new ArgumentException();
                        int N = InjectorsBase.GetLength(1);
                        var check = MultidimensionalArray.Create(N, N);
                        for (int j = 0; j < Jbase; j++) {
                            check.Clear();
                            check.AccEye(1.0);
                            check.Acc(-1.0, InjectorsBase.ExtractSubArrayShallow(j, -1, -1));

                            if (check.InfNorm() != 0.0)
                                throw new ArgumentException();

                        }
                    }

#endif
                    var iLPar = maxDgBasis.GridDat.iLogicalCells; // cells of *parent* grid

                    for (int j = 0; j < Jagg; j++) { // loop over aggregate cells
                        Debug.Assert(ArrayTools.ListEquals(Ag2Pt[j], C2F[j]));

                        int[] compCell = Ag2Pt[j];
                        int I = compCell.Length;


                        int iRoot = -1;
                        double maxSize = -1.0;
                        for(int i = 0; i < I; i++) {
                            double sz = iLPar.GetCellVolume(compCell[i]);
                            if (sz <= 0.0)
                                throw new ArithmeticException("found cell with non-positive volume.");
                            if(sz > maxSize) {
                                iRoot = i;
                                maxSize = sz;
                            }
                        }

                        //var gdat = ((maxDgBasis.GridDat) as Classic.GridData);
                        //double[] Sizes = compCell.Select(jPart => gdat.iLogicalCells.GetCellVolume(jPart)).ToArray();


                        Injectors_iLevel[j] = MultidimensionalArray.Create(I, Np, Np);
                        if (I > 1) {
                            // compute extrapolation
                            int[,] CellPairs = new int[I - 1, 2];
                            int cnt = 0;
                            for (int i = 0; i < I; i++) {
                                if (i != iRoot) {
                                    CellPairs[cnt, 0] = compCell[iRoot];
                                    CellPairs[cnt, 1] = compCell[i];
                                    cnt++;
                                }
                            }
                            var ExpolMtx = MultidimensionalArray.Create(I, Np, Np);
                            maxDgBasis.GetExtrapolationMatrices(CellPairs, ExpolMtx.ExtractSubArrayShallow(new int[] { 1, 0, 0 }, new int[] { I - 1, Np - 1, Np - 1 }));
                            for (int n = 0; n < Np; n++) {
                                ExpolMtx[0, n, n] = 1.0;
                            }

                            // Compute intermediate mass matrix
                            var MMtemp = MultidimensionalArray.Create(Np, Np);
                            MMtemp.Multiply(1.0, ExpolMtx, ExpolMtx, 0.0, "nm", "iln", "ilm");

                            // orthonormalize
                            //try {
                                MMtemp.SymmetricLDLInversion(ortho, null); // ortho is output, will be overwritten

                            //} catch (ArithmeticException ae) {
                            //    PlotScheisse(gdat.Grid, compCell);
                            //}
                            Injectors_iLevel[j].Multiply(1.0, ExpolMtx, ortho, 0.0, "inm", "ink", "km");
                        } else {
                            Injectors_iLevel[j].ExtractSubArrayShallow(0, -1, -1).AccEye(1.0);
                        }

                        // base level injector
                        var injBase = InjectorsBase.ExtractSubArrayShallow(compCell[0], -1, -1);
                        injBase.Set(Injectors_iLevel[j].ExtractSubArrayShallow(0, -1, -1));
                        Debug.Assert(InjectorsBaseReady[compCell[0]]);

                        for (int i = 1; i < I; i++) {
                            InjectorsBaseReady[compCell[i]] = false;
                        }



                    }

                    return Injectors_iLevel;
                }
            }


            /// <summary>
            /// computes the injector for multigrid level 2 and higher 
            /// </summary>
            /// <seealso cref="BuildInjector_Lv1"/>
            private static MultidimensionalArray[] BuildInjector_Lv2andup(
                Basis maxDgBasis, int Np,
                MultidimensionalArray InjectorsBase, bool[] InjectorsBaseReady,
                int Jagg, int[][] Ag2Pt_Fine, int[][] C2F) {
                using (new FuncTrace()) {

                    MultidimensionalArray[] Injectors_iLevel;
                    Injectors_iLevel = new MultidimensionalArray[Jagg];

                    for (int j = 0; j < Jagg; j++) { // loop over aggregate cells

                        // find cell pairs
                        int[] compCell = C2F[j];
                        int I = compCell.Length;

                        int[] BaseCells = new int[I];
                        for (int i = 0; i < I; i++) { // loop over parts
                            int jFine = compCell[i];
                            int[] jBaseS = Ag2Pt_Fine[jFine];
                            int II = jBaseS.Length;

                            BaseCells[i] = -1;
                            for (int ii = 0; ii < II; ii++) {
                                int jBase = jBaseS[ii];
                                if (InjectorsBaseReady[jBase]) {
                                    BaseCells[i] = jBase;
                                    break;
                                }
                            }
                            if (BaseCells[i] < 0)
                                throw new ApplicationException("Error in algorithm/data structure.");
                        }
                        int[,] CellPairs = new int[I - 1, 2];
                        for (int i = 0; i < I - 1; i++) {
                            CellPairs[i, 0] = BaseCells[0];
                            CellPairs[i, 1] = BaseCells[i + 1];
                        }

                        Injectors_iLevel[j] = MultidimensionalArray.Create(I, Np, Np);
                        if (I > 1) {
                            // ++++++++++++++++++++++++++++++++++
                            // 'normal' aggregation cell
                            // ++++++++++++++++++++++++++++++++++

                            // compute extrapolation (with respect to base grid)
                            var ExpolMtxBase = MultidimensionalArray.Create(I, Np, Np);
                            maxDgBasis.GetExtrapolationMatrices(CellPairs, ExpolMtxBase.ExtractSubArrayShallow(new int[] { 1, 0, 0 }, new int[] { I - 1, Np - 1, Np - 1 }));
                            for (int n = 0; n < Np; n++) {
                                ExpolMtxBase[0, n, n] = 1.0;
                            }

                            // compute extrapolation (with respect to finer level)
                            var ExpolMtx = MultidimensionalArray.Create(I, Np, Np);
                            var inv_injBase_i = MultidimensionalArray.Create(Np, Np);
                            for (int i = 0; i < I; i++) {
                                var ExpolMtxBase_i = ExpolMtxBase.ExtractSubArrayShallow(i, -1, -1);
                                var ExpolMtx_i = ExpolMtx.ExtractSubArrayShallow(i, -1, -1);
                                var injBase_i = InjectorsBase.ExtractSubArrayShallow(BaseCells[i], -1, -1);
                                Debug.Assert(InjectorsBaseReady[BaseCells[i]] == true);
                                Debug.Assert(injBase_i.InfNorm() > 0);

                                injBase_i.TriangularInvert(inv_injBase_i);
                                Debug.Assert(inv_injBase_i.InfNorm() > 0);

                                ExpolMtx_i.GEMM(1.0, inv_injBase_i, ExpolMtxBase_i, 0.0);
                                Debug.Assert(ExpolMtx_i.InfNorm() > 0);
                            }

                            // Compute intermediate mass matrix
                            var MMtemp = MultidimensionalArray.Create(Np, Np);
                            MMtemp.Multiply(1.0, ExpolMtx, ExpolMtx, 0.0, "nm", "iln", "ilm");

                            // orthonormalize
                            MultidimensionalArray ortho = MultidimensionalArray.Create(Np, Np);
                            MMtemp.SymmetricLDLInversion(ortho, null);
                            Injectors_iLevel[j].Multiply(1.0, ExpolMtx, ortho, 0.0, "inm", "ink", "km");
                        } else {
                            // +++++++++++++++++++++++++++++++++++++++++++++++++++++
                            // aggregate cell consists of only one cell 
                            // (degenerate case, may happen form time to time
                            // +++++++++++++++++++++++++++++++++++++++++++++++++++++
                            Injectors_iLevel[j].ExtractSubArrayShallow(0, -1, -1).AccEye(1.0);
                        }

                        // record injector to base grid
                        int jKeep = BaseCells[0];
                        var injBase_jKeep = InjectorsBase.ExtractSubArrayShallow(jKeep, -1, -1);
                        var injBase_jKeep_prev = injBase_jKeep.CloneAs();
                        var injLevel = Injectors_iLevel[j].ExtractSubArrayShallow(0, -1, -1);
                        injBase_jKeep.GEMM(1.0, injBase_jKeep_prev, injLevel, 0.0);

                        for (int i = 1; i < I; i++) {
                            int jDump = BaseCells[i];
                            Debug.Assert(InjectorsBaseReady[jDump] == true);
                            InjectorsBaseReady[jDump] = false;
#if DEBUG
                            InjectorsBase.ExtractSubArrayShallow(jDump, -1, -1).Clear();
#endif
                        }

                    }

                    return Injectors_iLevel;
                }
            }




            MultidimensionalArray CA(int _jAgg, int Np) {
                AggregationGridData ag = this.m_owner;
                var compCell = ag.iLogicalCells.AggregateCellToParts[_jAgg];
                int thisMgLevel = ag.MgLevel;


                var scl = m_owner.AncestorGrid.ChefBasis.Scaling;


                var R = MultidimensionalArray.Create(compCell.Length, Np, Np);
                for (int i = 0; i < compCell.Length; i++) {
                    int jG = compCell[i];
                    if (!m_owner.AncestorGrid.Cells.IsCellAffineLinear(jG))
                        throw new NotImplementedException("nonlin cell -- todo");
                    for (int n = 0; n < Np; n++) {
                        R[i, n, n] = scl[jG];
                    }
                }

#if DEBUG
                bool[] btouch = new bool[compCell.Length];
#endif

                int[] AggIndex = new int[] { _jAgg };
                _BasisData basisLevel = this;
                for (int mgLevelIdx = thisMgLevel; mgLevelIdx >= 0; mgLevelIdx--) {
                    AggregationGridData mgLevel = basisLevel.m_owner;
                    int[][] agg2part_parrent = mgLevel.ParentGrid.iLogicalCells.AggregateCellToParts;
#if DEBUG
                    btouch.Clear();
#endif

                    foreach (int jAgg in AggIndex) {
                        //MultidimensionalArray Inj_j;
                        //if (mgLevelIdx > 0) {
                        var Inj_j = basisLevel.Injectors[jAgg];
                        //} else {
                        //    Inj_j = MultidimensionalArray.Create(1, Np, Np);
                        //    for (int n = 0; n < Np; n++) {
                        //        m_CompositeBasis[jAgg][0, n, n] = 1.0;
                        //    }
                        //}


                        int[] FineAgg = mgLevel.jCellCoarse2jCellFine[jAgg];
                        Debug.Assert(FineAgg.Length == Inj_j.GetLength(0));

                        for (int iSrc = 0; iSrc < FineAgg.Length; iSrc++) { // loop over finer level cells
                            int jAgg_fine = FineAgg[iSrc];
                            // Inj_j[iSrc,-,-] is injector 
                            //   from cell 'jAgg' on level 'mgLevelIdx'      (coarse level)
                            //   to cell 'jAgg_fine' on level 'mgLevelIdx - 1' (fine level)

                            var Inj_j_iSrc = Inj_j.ExtractSubArrayShallow(iSrc, -1, -1);

                            int[] TargCells;// = mgLevel.ParentGrid.iLogicalCells.AggregateCellToParts[jAgg_fine];
                            if (agg2part_parrent != null)
                                TargCells = agg2part_parrent[jAgg_fine];
                            else
                                TargCells = new int[] { jAgg_fine };

                            foreach (int j in TargCells) {
                                int iTarg = Array.IndexOf(compCell, j);
                                if (iTarg < 0)
                                    throw new ApplicationException("error in alg");
#if DEBUG
                                if (btouch[iTarg] == true)
                                    throw new ApplicationException();
                                btouch[iTarg] = true;
#endif
                                var R_iTarg = R.ExtractSubArrayShallow(iTarg, -1, -1);

                                R_iTarg.Multiply(1.0, Inj_j_iSrc, R_iTarg.CloneAs(), 0.0, "nm", "nk", "km");

                                //if (thisMgLevel == 1 && Np == 10) {
                                //    var check = MultidimensionalArray.Create(Np, Np);
                                //    check.GEMM(1.0, R_iTarg, R_iTarg.Transpose(), 0.0);
                                //    check.AccEye(-1.0);
                                //    var bla = check.InfNorm();
                                //    Console.WriteLine("Check norm: " + bla);
                                //}
                            }
                        }
                    }


                    // Rekursions-Scheisse:
                    // - - - - - - - - - - -
                    //if (mgLevelIdx > 0) {
                    {
                        List<int> nextAggIndex = new List<int>();
                        foreach (int jAgg in AggIndex) {
                            int[] NextLevel = mgLevel.jCellCoarse2jCellFine[jAgg];
#if DEBUG
                            foreach (int i in NextLevel) {
                                Debug.Assert(nextAggIndex.Contains(i) == false);
                            }
#endif
                            nextAggIndex.AddRange(NextLevel);
                        }
                        AggIndex = nextAggIndex.ToArray();
                        if (m_owner.ParentGrid is AggregationGridData)
                            basisLevel = ((AggregationGridData)(m_owner.ParentGrid)).m_ChefBasis;
                        else
                            basisLevel = null;
                    }
                    //else {
                    //    AggIndex = null;
                    //    mgLevel = null;
                    //}
                }

                return R;
            }



        }



    }
}
