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
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;

namespace BoSSS.Solution.Multigrid {


    /// <summary>
    /// DG basis on an aggregation grid (<see cref="AggregationGrid"/>).
    /// </summary>
    public class AggregationGridBasis {


        /*
         * old version, inaccurate
         * 
         * 
        public static AggregationGridBasis[][] CreateSequence(IEnumerable<AggregationGrid> _agSeq, IEnumerable<Basis> dgBasisS) {

            // check input
            // -----------

            Basis maxDgBasis = null;
            XDGBasis maxXdgBasis = null;
            bool[] UseX = new bool[dgBasisS.Count()];
            bool AnyX = false, AnyNonX = false;
            for(int iBs = 0; iBs < dgBasisS.Count(); iBs++) {
                Basis b = dgBasisS.ElementAt(iBs);
                if(!object.ReferenceEquals(_agSeq.First().AncestorGrid, b.GridDat))
                    throw new ArgumentException("Mismatch between DG basis grid and multi-grid ancestor.");

                if(b is XDGBasis) {
                    var xb = b as XDGBasis;
                    if(maxDgBasis == null || b.Degree > maxDgBasis.Degree) {
                        maxDgBasis = xb.NonX_Basis;
                    }
                    if(maxXdgBasis == null || b.Degree > maxXdgBasis.Degree) {
                        maxXdgBasis = xb;
                    }
                    UseX[iBs] = true;
                    AnyX = true;
                } else {
                    if(maxDgBasis == null || b.Degree > maxDgBasis.Degree) {
                        maxDgBasis = b;
                    }
                    AnyNonX = true;
                }
            }

            AggregationGrid[] agSeq = _agSeq.ToArray();
            if (agSeq.Length <= 0)
                throw new ArgumentException();
            if (!object.ReferenceEquals(agSeq[0].ParentGrid, agSeq[0].AncestorGrid))
                throw new ArgumentException("Parent and Ancestor of 0th grid level must be equal.");
            GridData baseGrid = (GridData)(agSeq[0].AncestorGrid);
            
            for(int iLevel = 0; iLevel < agSeq.Length; iLevel++) {
                if (agSeq[iLevel].MgLevel != iLevel)
                    throw new ArgumentException("Grid levels must be provided in order.");
                if (!object.ReferenceEquals(agSeq[iLevel].AncestorGrid, baseGrid))
                    throw new ArgumentException("Mismatch in ancestor grid.");
                if(iLevel > 0) {
                    if(!object.ReferenceEquals(agSeq[iLevel].ParentGrid, agSeq[iLevel - 1])) {
                        throw new ArgumentException("Mismatch in parent grid at level " + iLevel + ".");
                    }
                }
            }

            if (baseGrid.CellPartitioning.TotalLength != agSeq[0].CellPartitioning.TotalLength)
                throw new ArgumentException("Mismatch in number of cells for level 0.");

            if (!object.ReferenceEquals(baseGrid, maxDgBasis.GridDat))
                throw new ArgumentException("Mismatch between DG basis grid and multi-grid ancestor.");

            // Project Bounding-Box basis
            // --------------------------
            var BB = baseGrid.GlobalBoundingBox;
            int Jbase = baseGrid.Cells.NoOfLocalUpdatedCells;
            int p = maxDgBasis.Degree;
            int Np = maxDgBasis.Length;
            PolynomialList polyList = maxDgBasis.Polynomials[0];

            MultidimensionalArray a = MultidimensionalArray.Create(Jbase, Np, Np);
            {
                CellQuadrature.GetQuadrature(new int[2] { Np, Np }, baseGrid,
                    (new CellQuadratureScheme()).Compile(baseGrid, p * 2),
                    delegate (int j0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                        NodeSet nodes = QR.Nodes;
                        int D = nodes.SpatialDimension;
                        NodeSet GlobalNodes = new NodeSet(baseGrid.Grid.RefElements[0], Length * nodes.NoOfNodes, D);
                        baseGrid.TransformLocal2Global(nodes, j0, Length, GlobalNodes.ResizeShallow(Length, nodes.NoOfNodes, D));

                        for (int d = 0; d < D; d++) {
                            var Cd = GlobalNodes.ExtractSubArrayShallow(-1, d);
                            Cd.Scale(2 / (BB.Max[d] - BB.Min[d]));
                            Cd.AccConstant((BB.Max[d] + BB.Min[d]) / (BB.Min[d] - BB.Max[d]));
                        }
#if DEBUG
                    Debug.Assert(GlobalNodes.Min() >= -1.00001);
                    Debug.Assert(GlobalNodes.Max() <= +1.00001);
#endif
                    GlobalNodes.LockForever();

                        var BasisValues = maxDgBasis.CellEval(nodes, j0, Length);
                        var PolyVals = MultidimensionalArray.Create(GlobalNodes.NoOfNodes, Np);
                        polyList.Evaluate(GlobalNodes, PolyVals);
                        PolyVals = PolyVals.ResizeShallow(Length, nodes.NoOfNodes, Np);

                        _EvalResult.Multiply(1.0, BasisValues, PolyVals, 0.0, "jknm", "jkn", "jkm");
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // _SaveIntegrationResults:
                    a.ExtractSubArrayShallow(new[] { i0, 0, 0 }, new[] { i0 + Length - 1, Np - 1, Np - 1 }).Acc(1.0, ResultsOfIntegration);
                    }).Execute();
            }

            // Compute intermediate mass matrix of bounding box basis
            // ------------------------------------------------------
            MultidimensionalArray Mass_Level = MultidimensionalArray.Create(Jbase, Np, Np);
            Mass_Level.Multiply(1.0, a, a, 0.0, "jlk", "jnl", "jnk");
            //a.ResizeShallow(Jbase * Np, Np).SaveToTextFile("c:\\tmp\\a.txt");
            //Mass_Level.ResizeShallow(Jbase * Np, Np).SaveToTextFile("c:\\tmp\\mass.txt");

            // Orthonormalize and create injection operator
            // --------------------------------------------

            MultidimensionalArray B_prevLevel, B_Level;

            // level 0
            int Jagg = agSeq[0].iLogicalCells.NoOfLocalUpdatedCells;
            B_Level = MultidimensionalArray.Create(Jagg, Np, Np);
            int[][] Ag2Pt = agSeq[0].iLogicalCells.AggregateCellToParts;
            for(int j = 0; j < Jagg; j++) {
                int jGeom;
                 if(Ag2Pt == null || Ag2Pt[j] == null) {
                    jGeom = j;
                } else {
                    if (Ag2Pt[j].Length != 1)
                        throw new ArgumentException();
                    jGeom = Ag2Pt[j][0];
                }
                if (jGeom != j)
                    throw new NotSupportedException("todo");

                a.ExtractSubArrayShallow(jGeom, -1, -1).InvertTo(B_Level.ExtractSubArrayShallow(j, -1, -1));
                //a.ExtractSubArrayShallow(jGeom, -1, -1).TriangularInvert(B_Level.ExtractSubArrayShallow(j, -1, -1));
            }

            MultidimensionalArray[][] Injectors = new MultidimensionalArray[agSeq.Length][];
            //var a_Level = a;

            // all other levels
            MultidimensionalArray Mass_prevLevel;
            for(int iLevel = 1; iLevel < agSeq.Length; iLevel++) { // loop over levels...
                B_prevLevel = B_Level;
                //a_prevLevel = a_Level;
                Mass_prevLevel = Mass_Level;
                Jagg = agSeq[iLevel].iLogicalCells.NoOfLocalUpdatedCells;
                Ag2Pt = agSeq[iLevel].iLogicalCells.AggregateCellToParts;
                int[][] C2F = agSeq[iLevel].jCellCoarse2jCellFine;

                var Injectors_iLevel = new MultidimensionalArray[Jagg];
                Injectors[iLevel] = Injectors_iLevel;

                B_Level = MultidimensionalArray.Create(Jagg, Np, Np);
                Mass_Level = MultidimensionalArray.Create(Jagg, Np, Np);
                for (int j = 0; j < Jagg; j++) { // loop over aggregate cells

                    var Mass_j = Mass_Level.ExtractSubArrayShallow(j, -1, -1);
                    foreach (int jF in C2F[j]) { // loop over aggregated cells
                        Mass_j.Acc(1.0, Mass_prevLevel.ExtractSubArrayShallow(jF, -1, -1));
                        //a_j.Acc(1.0, a_prevLevel.ExtractSubArrayShallow(jF, -1, -1));
                    }



                    // perform ortho-normalization:
                    var B_Level_j = B_Level.ExtractSubArrayShallow(j, -1, -1);
                    Mass_j.SymmetricLDLInversion(B_Level_j, default(double[]));
                    //a_j.TriangularInvert(B_Level_j, true);
#if DEBUG
                    {
                       

                        // assert that B_Level_j is upper-triangular
                        for (int n = 0; n < Np; n++) {
                            for (int m = n + 1; m < Np; m++) {
                                Debug.Assert(B_Level_j[m, n] == 0.0);
                            }
                        }
                    }
#endif
                    // invert 'B' from previous level
                    var invB_prevLevel = MultidimensionalArray.Create(B_prevLevel.Lengths);
                    for(int jP = 0; jP < B_prevLevel.GetLength(0); jP++) {
                        B_prevLevel.ExtractSubArrayShallow(jP, -1, -1).InvertTo(invB_prevLevel.ExtractSubArrayShallow(jP, -1, -1));
                    }
                    
                    // create injector
                    Injectors_iLevel[j] = MultidimensionalArray.Create(C2F[j].Length, Np, Np);
                    for (int l = 0; l < C2F[j].Length; l++) { // loop fine cells which make up the aggregate cell 
                        int jF = C2F[j][l];

                        var Inj_jl = Injectors_iLevel[j].ExtractSubArrayShallow(l, -1, -1);
                        var invB = invB_prevLevel.ExtractSubArrayShallow(jF, -1, -1);
                        Inj_jl.Multiply(1.0, invB, B_Level_j, 0.0, "nm", "nk", "km");
                    }
                }
            }

            // create basis sequence
            // ---------------------
            {
                AggregationGridBasis[][] ret = new AggregationGridBasis[agSeq.Length][];
                AggregationGridBasis[] Abs = new AggregationGridBasis[agSeq.Length];
                XdgAggregationBasis[] XAbs = new XdgAggregationBasis[agSeq.Length];
                for (int iLevel = 0; iLevel < agSeq.Length; iLevel++) {
                    if (AnyNonX)
                        Abs[iLevel] = new AggregationGridBasis(maxDgBasis, iLevel > 0 ? Abs[iLevel - 1] : null, agSeq[iLevel], Injectors[iLevel]);
                    if (AnyX)
                        XAbs[iLevel] = new XdgAggregationBasis(maxXdgBasis, iLevel > 0 ? XAbs[iLevel - 1] : null, agSeq[iLevel], Injectors[iLevel]);

                    ret[iLevel] = new AggregationGridBasis[UseX.Length];
                    for (int i = 0; i < UseX.Length; i++)
                        ret[iLevel][i] = UseX[i] ? XAbs[iLevel] : Abs[iLevel];
                }
                return ret;
            }
        }
        */


        /// <summary>
        /// Creation of a sequence of aggregation basis objects.
        /// </summary>
        /// <param name="_agSeq"></param>
        /// Sequence of aggregation grids.
        /// <param name="dgBasisS">
        /// List of DG basis objects on base grid.
        /// </param>
        /// <returns>
        /// - 1st index: correlates with grid level, i.e. with <paramref name="_agSeq"/>.
        /// - 2nd index: correlates with DG basis, <paramref name="dgBasisS"/>.
        /// </returns>
        public static AggregationGridBasis[][] CreateSequence(IEnumerable<AggregationGrid> _agSeq, IEnumerable<Basis> dgBasisS) {

            if (_agSeq.Count() <= 0)
                return new AggregationGridBasis[0][];
                    
            // check input
            // -----------

            Basis maxDgBasis = null;
            XDGBasis maxXdgBasis = null;
            bool[] UseX = new bool[dgBasisS.Count()];
            bool AnyX = false, AnyNonX = false;
            for(int iBs = 0; iBs < dgBasisS.Count(); iBs++) {
                Basis b = dgBasisS.ElementAt(iBs);
                if(!object.ReferenceEquals(_agSeq.First().AncestorGrid, b.GridDat))
                    throw new ArgumentException("Mismatch between DG basis grid and multi-grid ancestor.");

                if(b is XDGBasis) {
                    var xb = b as XDGBasis;
                    if(maxDgBasis == null || b.Degree > maxDgBasis.Degree) {
                        maxDgBasis = xb.NonX_Basis;
                    }
                    if(maxXdgBasis == null || b.Degree > maxXdgBasis.Degree) {
                        maxXdgBasis = xb;
                    }
                    UseX[iBs] = true;
                    AnyX = true;
                } else {
                    if(maxDgBasis == null || b.Degree > maxDgBasis.Degree) {
                        maxDgBasis = b;
                    }
                    AnyNonX = true;
                }
            }

            AggregationGrid[] agSeq = _agSeq.ToArray();
            if (agSeq.Length <= 0)
                throw new ArgumentException();
            if (!object.ReferenceEquals(agSeq[0].ParentGrid, agSeq[0].AncestorGrid))
                throw new ArgumentException("Parent and Ancestor of 0th grid level must be equal.");
            GridData baseGrid = (GridData)(agSeq[0].AncestorGrid);
            
            for(int iLevel = 0; iLevel < agSeq.Length; iLevel++) {
                if (agSeq[iLevel].MgLevel != iLevel)
                    throw new ArgumentException("Grid levels must be provided in order.");
                if (!object.ReferenceEquals(agSeq[iLevel].AncestorGrid, baseGrid))
                    throw new ArgumentException("Mismatch in ancestor grid.");
                if(iLevel > 0) {
                    if(!object.ReferenceEquals(agSeq[iLevel].ParentGrid, agSeq[iLevel - 1])) {
                        throw new ArgumentException("Mismatch in parent grid at level " + iLevel + ".");
                    }
                }
            }

            if (baseGrid.CellPartitioning.TotalLength != agSeq[0].CellPartitioning.TotalLength)
                throw new ArgumentException("Mismatch in number of cells for level 0.");

            if (!object.ReferenceEquals(baseGrid, maxDgBasis.GridDat))
                throw new ArgumentException("Mismatch between DG basis grid and multi-grid ancestor.");

            // Project Bounding-Box basis
            // --------------------------
            var BB = baseGrid.GlobalBoundingBox;
            int Jbase = baseGrid.Cells.NoOfLocalUpdatedCells;
            int p = maxDgBasis.Degree;
            int Np = maxDgBasis.Length;
            PolynomialList polyList = maxDgBasis.Polynomials[0];



            // -----------------------------------------------------------------------------------
            // injectors to upper level
            // - 1st index: level index; 0-th entry will always be null.
            // - 2nd index: aggregation cell index
            // For each multidimensional array:
            // - 1st index: enumeration of parts of upper level
            // - 2nd index: row index
            // - 3rd index: column index
            MultidimensionalArray[][] Injectors = new MultidimensionalArray[agSeq.Length][];
            // -----------------------------------------------------------------------------------

            MultidimensionalArray InjectorsBase = MultidimensionalArray.Create(Jbase, Np, Np);
            bool[] InjectorsBaseReady = new bool[Jbase];

            // level 0
            {
                int Jagg = agSeq[0].iLogicalCells.NoOfLocalUpdatedCells;
                int[][] Ag2Pt = agSeq[0].iLogicalCells.AggregateCellToParts;
                Debug.Assert(Jagg == Jbase);

                for (int j = 0; j < Jagg; j++) {
                    int jGeom;
                    if (Ag2Pt == null || Ag2Pt[j] == null) {
                        jGeom = j;
                    } else {
                        if (Ag2Pt[j].Length != 1)
                            throw new ArgumentException();
                        jGeom = Ag2Pt[j][0];
                    }
                    if (jGeom != j)
                        throw new NotSupportedException("todo");

                    InjectorsBase.ExtractSubArrayShallow(j, -1, -1).AccEye(1.0);
                    InjectorsBaseReady[j] = true;
                }
            }


            // level 1
            if(agSeq.Length >= 2) {
                int iLevel = 1;

                int Jagg = agSeq[iLevel].iLogicalCells.NoOfLocalUpdatedCells;
                int[][] Ag2Pt = agSeq[iLevel].iLogicalCells.AggregateCellToParts;
                int[][] C2F = agSeq[iLevel].jCellCoarse2jCellFine;

                Debug.Assert(Ag2Pt.Length >= Jagg);
                Debug.Assert(C2F.Length >= Jagg);

#if DEBUG
                var InjectorsBase_clone = InjectorsBase.CloneAs();
                var InjectorsBaseReady_clone = InjectorsBaseReady.CloneAs();
                int[][] Ag2Pt_Fine = Jbase.ForLoop(j => new int[1] { j });
                MultidimensionalArray[] InjCheck = BuildInjector_Lv2andup(maxDgBasis, Np, InjectorsBase_clone, InjectorsBaseReady_clone, Jagg, Ag2Pt_Fine, C2F);
#endif
                Injectors[iLevel] = BuildInjector_Lv1(maxDgBasis, Np, InjectorsBase, InjectorsBaseReady, Jagg, Ag2Pt, C2F);
#if DEBUG
                for(int jAgg = 0; jAgg < Math.Max(InjCheck.Length, Injectors[iLevel].Length); jAgg++) {
                    double err = InjCheck[jAgg].L2Dist( Injectors[iLevel][jAgg]);
                    Debug.Assert(err < Math.Max(InjCheck[jAgg].L2Norm(), Injectors[iLevel][jAgg].L2Norm()) * 1e-8);
                }
#endif
            } else {
                //ortho_Level = null;
            }


            // all other levels
            for (int iLevel = 2; iLevel < agSeq.Length; iLevel++) { // loop over levels...
                int Jagg = agSeq[iLevel].iLogicalCells.NoOfLocalUpdatedCells;
                int[][] Ag2Pt = agSeq[iLevel].iLogicalCells.AggregateCellToParts;
                int[][] Ag2Pt_Fine = agSeq[iLevel - 1].iLogicalCells.AggregateCellToParts;
                int[][] C2F = agSeq[iLevel].jCellCoarse2jCellFine;

                Debug.Assert(Ag2Pt.Length >= Jagg);
                Debug.Assert(C2F.Length >= Jagg);

               

                //ortho_PrevLevel = ortho_Level;
                //ortho_Level = MultidimensionalArray.Create(Jagg, Np, Np);

                Injectors[iLevel] = BuildInjector_Lv2andup(maxDgBasis, Np, InjectorsBase, InjectorsBaseReady, Jagg, Ag2Pt_Fine, C2F);
            }

            // create basis sequence
            // ---------------------
            {
                AggregationGridBasis[][] ret = new AggregationGridBasis[agSeq.Length][];
                AggregationGridBasis[] Abs = new AggregationGridBasis[agSeq.Length];
                XdgAggregationBasis[] XAbs = new XdgAggregationBasis[agSeq.Length];
                for (int iLevel = 0; iLevel < agSeq.Length; iLevel++) {
                    if (AnyNonX) {
                        Abs[iLevel] = new AggregationGridBasis(maxDgBasis, iLevel > 0 ? Abs[iLevel - 1] : null, agSeq[iLevel], Injectors[iLevel]);
                    }
                    if (AnyX) {
                        XAbs[iLevel] = new XdgAggregationBasis(maxXdgBasis, iLevel > 0 ? XAbs[iLevel - 1] : null, agSeq[iLevel], Injectors[iLevel]);
                    }

                    ret[iLevel] = new AggregationGridBasis[UseX.Length];
                    for (int i = 0; i < UseX.Length; i++)
                        ret[iLevel][i] = UseX[i] ? XAbs[iLevel] : Abs[iLevel];
                }
                return ret;
            }
        }

        /// <summary>
        /// computes the injector for multigrid level 2 and higher 
        /// </summary>
        /// <seealso cref="BuildInjector_Lv2andup"/>
        private static MultidimensionalArray[] BuildInjector_Lv1(
            Basis maxDgBasis, int Np, 
            MultidimensionalArray InjectorsBase, bool[] InjectorsBaseReady, 
            int Jagg, int[][] Ag2Pt, int[][] C2F) {
            using(new FuncTrace()) {
                MultidimensionalArray ortho = MultidimensionalArray.Create(Np, Np);

                MultidimensionalArray[] Injectors_iLevel = new MultidimensionalArray[Jagg];

                for(int j = 0; j < Jagg; j++) { // loop over aggregate cells
                    Debug.Assert(ArrayTools.ListEquals(Ag2Pt[j], C2F[j]));

                    int[] compCell = Ag2Pt[j];

                    int I = compCell.Length;

                    Injectors_iLevel[j] = MultidimensionalArray.Create(I, Np, Np);
                    if(I > 1) {
                        // compute extrapolation
                        int[,] CellPairs = new int[I - 1, 2];
                        for(int i = 0; i < I - 1; i++) {
                            CellPairs[i, 0] = compCell[0];
                            CellPairs[i, 1] = compCell[i + 1];
                        }
                        var ExpolMtx = MultidimensionalArray.Create(I, Np, Np);
                        maxDgBasis.GetExtrapolationMatrices(CellPairs, ExpolMtx.ExtractSubArrayShallow(new int[] { 1, 0, 0 }, new int[] { I - 1, Np - 1, Np - 1 }));
                        for(int n = 0; n < Np; n++) {
                            ExpolMtx[0, n, n] = 1.0;
                        }

                        // Compute intermediate mass matrix
                        var MMtemp = MultidimensionalArray.Create(Np, Np);
                        MMtemp.Multiply(1.0, ExpolMtx, ExpolMtx, 0.0, "nm", "iln", "ilm");

                        // orthonormalize
                        MMtemp.SymmetricLDLInversion(ortho, null);
                        Injectors_iLevel[j].Multiply(1.0, ExpolMtx, ortho, 0.0, "inm", "ink", "km");
                    } else {
                        Injectors_iLevel[j].ExtractSubArrayShallow(0, -1, -1).AccEye(1.0);
                    }

                    // base level injector
                    var injBase = InjectorsBase.ExtractSubArrayShallow(compCell[0], -1, -1);
                    injBase.Set(Injectors_iLevel[j].ExtractSubArrayShallow(0, -1, -1));
                    Debug.Assert(InjectorsBaseReady[compCell[0]]);

                    for(int i = 1; i < I; i++) {
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
            using(new FuncTrace()) {

                MultidimensionalArray[] Injectors_iLevel;
                Injectors_iLevel = new MultidimensionalArray[Jagg];

                for(int j = 0; j < Jagg; j++) { // loop over aggregate cells

                    // find cell pairs
                    int[] compCell = C2F[j];
                    int I = compCell.Length;

                    int[] BaseCells = new int[I];
                    for(int i = 0; i < I; i++) { // loop over parts
                        int jFine = compCell[i];
                        int[] jBaseS = Ag2Pt_Fine[jFine];
                        int II = jBaseS.Length;

                        BaseCells[i] = -1;
                        for(int ii = 0; ii < II; ii++) {
                            int jBase = jBaseS[ii];
                            if(InjectorsBaseReady[jBase]) {
                                BaseCells[i] = jBase;
                                break;
                            }
                        }
                        if(BaseCells[i] < 0)
                            throw new ApplicationException("Error in algorithm/data structure.");
                    }
                    int[,] CellPairs = new int[I - 1, 2];
                    for(int i = 0; i < I - 1; i++) {
                        CellPairs[i, 0] = BaseCells[0];
                        CellPairs[i, 1] = BaseCells[i + 1];
                    }

                    Injectors_iLevel[j] = MultidimensionalArray.Create(I, Np, Np);
                    if(I > 1) {
                        // ++++++++++++++++++++++++++++++++++
                        // 'normal' aggregation cell
                        // ++++++++++++++++++++++++++++++++++

                        // compute extrapolation (with respect to base grid)
                        var ExpolMtxBase = MultidimensionalArray.Create(I, Np, Np);
                        maxDgBasis.GetExtrapolationMatrices(CellPairs, ExpolMtxBase.ExtractSubArrayShallow(new int[] { 1, 0, 0 }, new int[] { I - 1, Np - 1, Np - 1 }));
                        for(int n = 0; n < Np; n++) {
                            ExpolMtxBase[0, n, n] = 1.0;
                        }

                        // compute extrapolation (with respect to finer level)
                        var ExpolMtx = MultidimensionalArray.Create(I, Np, Np);
                        var inv_injBase_i = MultidimensionalArray.Create(Np, Np);
                        for(int i = 0; i < I; i++) {
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

                    for(int i = 1; i < I; i++) {
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



        static GridData GetGridData(AggregationGrid ag) {
            if(ag.ParentGrid is GridData) {
                return ((GridData)ag.ParentGrid);
            } else if(ag.ParentGrid is AggregationGrid) {
                return GetGridData((AggregationGrid)ag.ParentGrid);
            } else {
                throw new NotSupportedException();
            }
        }


        public AggregationGridBasis ParentBasis {
            get;
            private set;
        }


        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="b">DG basis on original grid</param>
        /// <param name="ag"></param>
        /// <param name="Injection">injection operator</param>
        protected AggregationGridBasis(Basis b, AggregationGridBasis parentBasis, AggregationGrid ag, MultidimensionalArray[] Injection) {
            using (new FuncTrace()) {
                if (!object.ReferenceEquals(b.GridDat, GetGridData(ag)))
                    throw new ArgumentException("mismatch in grid data object.");
                this.DGBasis = b;
                this.AggGrid = ag;

                if (parentBasis != null) {
                    if (!object.ReferenceEquals(ag.ParentGrid, parentBasis.AggGrid))
                        throw new ArgumentException("mismatch in parent grid.");
                } else {
                    if (ag.MgLevel != 0)
                        throw new ArgumentException();

                }

                ParentBasis = parentBasis;
#if DEBUG
                if (ag.MgLevel != 0) {
                    if (Injection.Length != ag.iLogicalCells.NoOfLocalUpdatedCells)
                        throw new ArgumentException();
                }
#endif
                m_InjectionOperator = Injection;

                //SetupProlongationOperator();
            }
        }

        MultidimensionalArray[] m_InjectionOperator;


        /// <summary>
        /// Injection/prolongation operator to finer grid level
        /// - array index: aggregate cell index; 
        /// - 1st index into <see cref="MultidimensionalArray"/>: index within aggregation basis, correlates with 2nd index into <see cref="AggregationGrid.jCellCoarse2jCellFine"/>
        /// - 2nd index into <see cref="MultidimensionalArray"/>: row
        /// - 3rd index into <see cref="MultidimensionalArray"/>: column
        /// - content: local cell index into the original grid, see <see cref="Foundation.Grid.ILogicalCellData.AggregateCellToParts"/>
        /// </summary>
        public MultidimensionalArray[] InjectionOperator {
            get {
                return m_InjectionOperator;
            }
        }

<<<<<<< HEAD
                        MassMatrix.AccEye(-1.0);
                        Debug.Assert(MassMatrix.InfNorm() < 1.0e-8);
#endif
                    }
=======
>>>>>>> 194f8c566ec46e985d0a6f0b37e02f3ee621ccaa
        

        /// <summary>
        /// restricts/projects a vector from the full grid (<see cref="BoSSS.Foundation.Grid.Classic.GridData"/>)
        /// to the aggregated grid (<see cref="AggGrid"/>).
        /// </summary>
        /// <param name="FullGridVector">input;</param>
        /// <param name="AggGridVector">output;</param>
        virtual public void RestictFromFullGrid<T, V>(T FullGridVector, V AggGridVector)
            where T : IList<double>
            where V : IList<double> //
        {
            var fullMapping = new UnsetteledCoordinateMapping(this.DGBasis);
            if (FullGridVector.Count != fullMapping.LocalLength)
                throw new ArgumentException("mismatch in vector length", "FullGridVector");
            int L = this.LocalDim;
            if (AggGridVector.Count != L)
                throw new ArgumentException("mismatch in vector length", "AggGridVector");

            var ag = this.AggGrid;
            var agCls = ag.iLogicalCells.AggregateCellToParts;
            int JAGG = ag.iLogicalCells.NoOfLocalUpdatedCells;
            int N = this.DGBasis.Length;

            var Buffer = MultidimensionalArray.Create(agCls.Max(cc => cc.Length), N);
            var Aggcoords = MultidimensionalArray.Create(N);

            for (int jAgg = 0; jAgg < JAGG; jAgg++) { // loop over all aggregated cells...
                int[] agCl = agCls[jAgg];
                int K = agCl.Length;

                MultidimensionalArray coords = (K * N == Buffer.GetLength(0)) ? Buffer : Buffer.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { K - 1, N - 1 });

                for (int k = 0; k < K; k++) { // loop over the cells which form the aggregated cell...
                    int jCell = agCl[k];
                    int j0 = fullMapping.LocalUniqueCoordinateIndex(0, jCell, 0);
                    for (int n = 0; n < N; n++) {
                        coords[k, n] = FullGridVector[n + j0];
                    }
                }

                Aggcoords.Clear();
                Aggcoords.Multiply(1.0, this.CompositeBasis[jAgg], coords, 0.0, "n", "kmn", "km");

                int i0 = jAgg * N;
                for (int n = 0; n < N; n++)
                    AggGridVector[n + i0] = Aggcoords[n];
            }
        }


        /// <summary>
        /// Prolongates/injects a vector from the full grid (<see cref="BoSSS.Foundation.Grid.GridData"/>)
        /// to the aggregated grid (<see cref="AggGrid"/>).
        /// </summary>
        /// <param name="FullGridVector">output;</param>
        /// <param name="AggGridVector">input;</param>
        virtual public void ProlongateToFullGrid<T, V>(T FullGridVector, V AggGridVector)
            where T : IList<double>
            where V : IList<double> //
        {
            var fullMapping = new UnsetteledCoordinateMapping(this.DGBasis);
            if(FullGridVector.Count != fullMapping.LocalLength)
                throw new ArgumentException("mismatch in vector length", "FullGridVector");
            int L = this.LocalDim;
            if(AggGridVector.Count != L)
                throw new ArgumentException("mismatch in vector length", "AggGridVector");


            var ag = this.AggGrid;
            var agCls = ag.iLogicalCells.AggregateCellToParts;
            int JAGG = ag.iLogicalCells.NoOfLocalUpdatedCells;
            int N = this.DGBasis.Length;

            var FulCoords = new double[N];
            var AggCoords = new double[N];

            //MultidimensionalArray Trf = MultidimensionalArray.Create(N, N);

            for(int jAgg = 0; jAgg < JAGG; jAgg++) {
                int[] agCl = agCls[jAgg];
                int K = agCl.Length;

                int i0 = jAgg * N;
                for(int n = 0; n < N; n++)
                    AggCoords[n] = AggGridVector[n + i0];

                for(int k = 0; k < K; k++) { // loop over the cells wich form the aggregated cell...
                    int jCell = agCl[k];
                    int j0 = fullMapping.LocalUniqueCoordinateIndex(0, jCell, 0);

                    //Trf.Clear();
                    //Trf.Acc(1.0, this.CompositeBasis[jAgg].ExtractSubArrayShallow(k, -1, -1));
                    var Trf = this.CompositeBasis[jAgg].ExtractSubArrayShallow(k, -1, -1);

                    //Trf.Solve(FulCoords, AggCoords);
                    Trf.gemv(1.0, AggCoords, 0.0, FulCoords);

                    for(int n = 0; n < N; n++) {
                        FullGridVector[j0 + n] = FulCoords[n];
                    }
                }
            }
        }

        /// <summary>
        /// Restriction operator from the base grid to some multigrid level,
        /// for a single variable/DG field of a <see cref="MultigridMapping"/>.
        /// </summary>
        /// <param name="rest">Output</param>
        /// <param name="mgMap"></param>
        /// <param name="iFld">DG field index within <see cref="mgMap"/>.</param>
        /// <remarks>
        /// Not intended for direct user interaction, mainly used by
        /// used by <see cref="MultigridMapping.FromOtherLevelMatrix(MultigridMapping)"/>
        /// </remarks>
        virtual public void GetRestrictionMatrix(BlockMsrMatrix rest, MultigridMapping mgMap, int iFld) {
            if(!object.ReferenceEquals(mgMap.AggBasis[iFld], this))
                throw new ArgumentException();

            UnsetteledCoordinateMapping fullmap = mgMap.ProblemMapping;
            int[] degree = mgMap.DgDegree;

            foreach(var b in fullmap.BasisS) {
                if(!b.IsSubBasis(this.DGBasis))
                    throw new ArgumentException("");
                if(b.MaximalLength != b.MinimalLength)
                    throw new NotSupportedException();
            }
            if(fullmap.BasisS.Count != degree.Length)
                throw new ArgumentException();
            int NoFld = degree.Length;

            int[] FullLength = fullmap.BasisS.Select(b => b.Length).ToArray();
            int[] RestLength = NoFld.ForLoop(
                l => fullmap.BasisS[l].Polynomials[0].Where(poly => poly.AbsoluteDegree <= degree[l]).Count());
            int[] RestOffset = new int[NoFld];
            for(int f = 1; f < NoFld; f++) {
                RestOffset[f] = RestOffset[f - 1] + RestLength[f - 1];
            }

            int NR = RestLength.Sum();
            int JAGG = this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            var ag = this.AggGrid;
            var agCls = ag.iLogicalCells.AggregateCellToParts;


            //Partitioning RowMap = new Partitioning(JAGG * NR);
            int MpiOffset_row = rest.RowPartitioning.i0;
            
            for(int jagg = 0; jagg < JAGG; jagg++) {
                int[] agCl = agCls[jagg];
                int K = agCl.Length;

                for (int k = 0; k < K; k++) { // loop over the cells which form the aggregated cell...
                    int jCell = agCl[k];

                    int M = RestLength[iFld];
                    int N = FullLength[iFld];
                    var Block = this.CompositeBasis[jagg].ExtractSubArrayShallow(new int[] { k, 0, 0 }, new int[] { k - 1, N - 1, M - 1 });

                    //for(int f = 0; f < NoFld; f++) { // loop over DG fields in mapping...

                    int i0Col = fullmap.GlobalUniqueCoordinateIndex(iFld, jCell, 0);
                    int i0Row = jagg * NR + RestOffset[iFld] + MpiOffset_row;

                    //rest.AccBlock(i0Row, i0Col, 1.0, Block);
                    for (int m = 0; m < M; m++) {
                        for (int n = 0; n < N; n++) {
                            rest[i0Row + m, i0Col + n] = Block[n, m];
                        }
                    }

                    //if (ag.MgLevel == 1 && N == 10) {
                    //    //var check = MultidimensionalArray.Create(N, N);
                    //    var RRt = Block.GEMM(Block.Transpose());
                    //    var check = RRt.GEMM(RRt);

                    //    check.AccEye(-1.0);
                    //    var bla = check.InfNorm();
                    //    Console.WriteLine("Check norm: " + bla);
                    //}

                }
            }

            if (ag.MgLevel == 1 ) {
                
                var result = BlockMsrMatrix.Multiply(rest, rest.Transpose());
                
                var resultT = result.Transpose();
                BlockMsrMatrix ShoudBeId;
                if(result.RowPartitioning.TotalLength < result.ColPartition.TotalLength)
                    ShoudBeId = BlockMsrMatrix.Multiply(result, resultT);
                else
                    ShoudBeId = BlockMsrMatrix.Multiply(resultT, result);

                ShoudBeId.AccEyeSp(-1.0);

                double ShouldBeID_Norm = ShoudBeId.InfNorm();
                Debug.Assert(ShouldBeID_Norm < 1.0e-8);
                //Console.WriteLine("Id norm {0} \t (lävel {1})", ShouldBeID_Norm, this.AggGrid.MgLevel);

            }
        }

        
        public AggregationGrid AggGrid {
            get;
            private set;
        }

        public Basis DGBasis {
            get;
            private set;
        }

        public virtual int MaximalLength {
            get {
                return this.DGBasis.MaximalLength;
            }
        }

        public virtual int MinimalLength {
            get {
                return this.DGBasis.MinimalLength;
            }
        }


        int[] m_Lengths;

        public virtual int GetMaximalLength(int p) {
            Debug.Assert(this.DGBasis.MaximalLength == this.DGBasis.MinimalLength);
            return this.GetLength(0, p);
        }

        public virtual int GetMinimalLength(int p) {
            Debug.Assert(this.DGBasis.MaximalLength == this.DGBasis.MinimalLength);
            return this.GetLength(0, p);
        }
        
        public virtual int GetLength(int jCell, int p) {
            GetNp();

            return m_Lengths[p];
        }

       
        /// <summary>
        /// The projector in the L2 Norm, from the space defined by the basis <see cref="DGBasis"/> on the original,
        /// onto the DG space on the aggregate grid.
        /// - array index: aggregate cell index; 
        /// - 1st index into <see cref="MultidimensionalArray"/>: index within aggregation basis 
        /// - 2nd index into <see cref="MultidimensionalArray"/>: row
        /// - 3rd index into <see cref="MultidimensionalArray"/>: column
        /// - content: local cell index into the original grid, see <see cref="Foundation.Grid.ILogicalCellData.AggregateCellToParts"/>
        /// </summary>
        /// <remarks>
        /// This method does not scale linear with problem size, its only here for reference/testing purpose.
        /// </remarks>
        public MultidimensionalArray[] CompositeBasis {
            get {
                if(m_CompositeBasis == null) {
                    SetupCompositeBasis();
                }
                return m_CompositeBasis;
            }
        }
        
        MultidimensionalArray CA(int _jAgg) {
            AggregationGrid ag = this.AggGrid;
            var compCell = ag.iLogicalCells.AggregateCellToParts[_jAgg];
            int thisMgLevel = ag.MgLevel;
            int Np = this.DGBasis.Length;

            var R = MultidimensionalArray.Create(compCell.Length, Np, Np);
            for(int i = 0; i < compCell.Length; i++) {
                for(int n = 0; n < Np; n++) {
                    R[i, n, n] = 1.0;
                }
            }

#if DEBUG
            bool[] btouch = new bool[compCell.Length];
#endif

            int[] AggIndex = new int[] { _jAgg };
            AggregationGridBasis basisLevel = this;
            for (int mgLevelIdx = thisMgLevel; mgLevelIdx > 0; mgLevelIdx--) {
                AggregationGrid mgLevel = basisLevel.AggGrid;
#if DEBUG
                btouch.Clear();
#endif

                foreach (int jAgg in AggIndex) {
                    //MultidimensionalArray Inj_j;
                    //if (mgLevelIdx > 0) {
                    var Inj_j = basisLevel.InjectionOperator[jAgg];
                    //} else {
                    //    Inj_j = MultidimensionalArray.Create(1, Np, Np);
                    //    for (int n = 0; n < Np; n++) {
                    //        m_CompositeBasis[jAgg][0, n, n] = 1.0;
                    //    }
                    //}


                    int[] FineAgg = mgLevel.jCellCoarse2jCellFine[jAgg];
                    Debug.Assert(FineAgg.Length == Inj_j.GetLength(0));

                    for(int iSrc = 0; iSrc < FineAgg.Length; iSrc++) { // loop over finer level cells
                        int jAgg_fine = FineAgg[iSrc];
                        // Inj_j[iSrc,-,-] is injector 
                        //   from cell 'jAgg' on level 'mgLevelIdx'      (coarse level)
                        //   to cell 'jAgg_fine' on level 'mgLevelIdx - 1' (fine level)

                        var Inj_j_iSrc = Inj_j.ExtractSubArrayShallow(iSrc, -1, -1);

                        int[] TargCells = mgLevel.ParentGrid.iLogicalCells.AggregateCellToParts[jAgg_fine];
                        
                        foreach (int j in TargCells) {
                            int iTarg = Array.IndexOf(compCell, j);
                            if (iTarg < 0)
                                throw new ApplicationException("error in alg");
#if DEBUG
                            if(btouch[iTarg] == true)
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
                    foreach(int jAgg in AggIndex) {
                        int[] NextLevel = mgLevel.jCellCoarse2jCellFine[jAgg];
#if DEBUG
                        foreach(int i in NextLevel) {
                            Debug.Assert(nextAggIndex.Contains(i) == false);
                        }
#endif
                        nextAggIndex.AddRange(NextLevel);
                    }
                    AggIndex = nextAggIndex.ToArray();
                    basisLevel = basisLevel.ParentBasis;
                } 
                //else {
                //    AggIndex = null;
                //    mgLevel = null;
                //}
            }

            return R;
        }


        MultidimensionalArray[] m_CompositeBasis;

        void SetupCompositeBasis() {
            using(new FuncTrace()) {
                Basis b = this.DGBasis;
                AggregationGrid ag = this.AggGrid;
                Debug.Assert(object.ReferenceEquals(b.GridDat, GetGridData(ag)));
                int N = b.Length;
                

                int JAGG = ag.iLogicalCells.NoOfLocalUpdatedCells;
                m_CompositeBasis = new MultidimensionalArray[JAGG];

                for(int jAgg = 0; jAgg < JAGG; jAgg++) { // loop over agglomerated cells...

                    m_CompositeBasis[jAgg] = CA(jAgg);

                    /*
                    var compCell = ag.iLogicalCells.AggregateCellToParts[jAgg];
                    
                    if(compCell.Length == 1) {


                        m_CompositeBasis[jAgg] = MultidimensionalArray.Create(1, N, N);
                        for(int n = 0; n < N; n++) {
                            m_CompositeBasis[jAgg][0, n, n] = 1.0;
                        }

                    } else {
                        // compute extrapolation basis
                        // ===========================

                        int I = compCell.Length - 1;
                        int[,] CellPairs = new int[I, 2];

                        for(int i = 0; i < I; i++) {
                            CellPairs[i, 0] = compCell[0];
                            CellPairs[i, 1] = compCell[i + 1];
                        }
                        var ExpolMtx = MultidimensionalArray.Create(I + 1, N, N);
                        b.GetExtrapolationMatrices(CellPairs, ExpolMtx.ExtractSubArrayShallow(new int[] { 1, 0, 0 }, new int[] { I, N - 1, N - 1 }));
                        for(int n = 0; n < N; n++) {
                            ExpolMtx[0, n, n] = 1.0;
                        }

                        // compute mass matrix
                        // ===================

                        var MassMatrix = MultidimensionalArray.Create(N, N); // intermediate mass matrix
                        MassMatrix.Multiply(1.0, ExpolMtx, ExpolMtx, 0.0, "lm", "kim", "kil");


                        // change to orthonormal basis
                        // ===========================
                        MultidimensionalArray B = MultidimensionalArray.Create(N, N);
                        MassMatrix.SymmetricLDLInversion(B, default(double[]));


                        m_CompositeBasis[jAgg] = MultidimensionalArray.Create(ExpolMtx.Lengths);
                        m_CompositeBasis[jAgg].Multiply(1.0, ExpolMtx, B, 0.0, "imn", "imk", "kn");

                        // check
                        // =====
#if DEBUG
                        MassMatrix.Clear();
                        for(int k = 0; k <= I; k++) {
                            for(int l = 0; l < N; l++) { // over rows of mass matrix ...
                                for(int m = 0; m < N; m++) { // over columns of mass matrix ...

                                    double mass_lm = 0.0;

                                    for(int i = 0; i < N; i++) {
                                        mass_lm += m_CompositeBasis[jAgg][k, i, m] * m_CompositeBasis[jAgg][k, i, l];
                                    }

                                    MassMatrix[l, m] += mass_lm;
                                }
                            }
                        }

                        MassMatrix.AccEye(-1.0);
                        Debug.Assert(MassMatrix.InfNorm() < 1.0e-9);
#endif
                    }
                    //*/
                }
            }
        }



        /// <summary>
        /// Returns a mapping form polynomial degree to DOF per cell.
        /// </summary>
        /// <returns></returns>
        public int[] GetNp() {
            if(m_Lengths == null) {
                m_Lengths = new int[this.DGBasis.Degree + 1];
                for(int pp = 0; pp < m_Lengths.Length; pp++) {
                    m_Lengths[pp] = this.DGBasis.Polynomials[0].Where(poly => poly.AbsoluteDegree <= pp).Count();
                }
            }
            return m_Lengths.CloneAs();
        }


        /// <summary>
        /// Local vector-space dimension.
        /// </summary>
        virtual public int LocalDim {
            get {
                return this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells * this.DGBasis.Length;
            }
        }


        /// <summary>
        /// for XDG, the cell mode index <paramref name="n"/> may not be equal
        /// in the full and the aggregated grid. This method performs the transformation.
        /// </summary>
        virtual internal int N_Murks(int j, int n, int N) {
            Debug.Assert(j >= 0 && j < this.DGBasis.GridDat.iLogicalCells.Count);
            Debug.Assert(n >= 0 && n < this.DGBasis.GetLength(j));
            Debug.Assert(n < N);
            Debug.Assert(this.AggGrid.iLogicalCells.NoOfLocalUpdatedCells == this.DGBasis.GridDat.iLogicalCells.Count);
            return n;
        }

        virtual internal bool ReqModeIndexTrafo {
            get {
                return false;
            }
        }

        int[][] m_ModeIndexForDegree;

        virtual public int[] ModeIndexForDegree(int j, int p, int Pmax) {
            if(m_ModeIndexForDegree == null) {
                int Psupermax = this.DGBasis.Degree;
                m_ModeIndexForDegree = new int[Psupermax + 1][];
                for(int _p = 0; _p <= Psupermax; _p++) {
                    m_ModeIndexForDegree[_p] = new int[0];
                }

                Debug.Assert(this.DGBasis.Polynomials.Count == 1);
                var polys = this.DGBasis.Polynomials[0];

                for(int n = 0; n < polys.Count; n++) {
                    n.AddToArray(ref m_ModeIndexForDegree[polys[n].AbsoluteDegree]);
                }
            }
            return m_ModeIndexForDegree[p];
        }
        
    }

}
