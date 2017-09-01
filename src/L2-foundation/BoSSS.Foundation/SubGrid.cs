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
using System.Collections;
using System.Collections.Generic;
using BoSSS.Platform;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;

namespace BoSSS.Foundation.Grid {

    /// <summary>
    /// denoted a subset of the computational grid
    /// </summary>
    public class SubGrid {

        IGridData m_GridData;

        /// <summary>
        /// the grid that this object is associated with
        /// </summary>
        public IGridData GridData {
            get {
                return m_GridData;
            }
        }

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="volMask">
        /// volume mask for those cells which should be contained in the subgrid
        /// </param>
        public SubGrid(CellMask volMask) {
            MPICollectiveWatchDog.Watch();
            this.m_VolumeMask = volMask;
            m_GridData = volMask.GridData;
        }

        private int[] m_LocalCellIndex2SubgridIndex;

        /// <summary>
        /// mapping from local (within the current MPI process) cell indices to
        /// the subgrid;<br/>
        /// index: a local cell index <em>j</em>, i.e. 0 &#8804; <em>j</em> &lt;
        /// <see cref="GridData.CellData.NoOfCells"/>;<br/>
        /// content: a subgrid index <em>k</em>, i.e. 0 &#8804; <em>k</em> &lt;
        /// <see cref="LocalNoOfCells"/>, or a negative number if cell
        /// <em>j</em> is not contained in the subgrid.;
        /// </summary>
        /// <remarks>
        /// The inverse of this mapping is
        /// <see cref="SubgridIndex2LocalCellIndex"/>;<br/>
        /// For performance reasons, this array is only computed once, and
        /// never cloned;
        /// </remarks>
        public int[] LocalCellIndex2SubgridIndex {
            get {
                if (m_LocalCellIndex2SubgridIndex == null) {
                    m_LocalCellIndex2SubgridIndex = ComputeLocalCellIndex2SubgridIndex();
                }
                return m_LocalCellIndex2SubgridIndex;
            }
        }

        /// <summary>
        /// Computes <see cref="LocalCellIndex2SubgridIndex"/>. Invokes
        /// <see cref="ComputeSubgridIndex2LocalCellIndex"/> which requires MPI
        /// communication.
        /// </summary>
        /// <returns></returns>
        private int[] ComputeLocalCellIndex2SubgridIndex() {
            int[] inverse = SubgridIndex2LocalCellIndex;
            int JE = m_GridData.iLogicalCells.NoOfCells;
            int[] localCellIndex2SubgridIndex = new int[JE];
            ArrayTools.SetAll(localCellIndex2SubgridIndex, int.MinValue);

            int K = inverse.Length;
            for (int k = 0; k < K; k++) {
                localCellIndex2SubgridIndex[inverse[k]] = k;
            }

            return localCellIndex2SubgridIndex;
        }

        private int[] m_SubgridIndex2LocalCellIndex;

        /// <summary>
        /// mapping from subgrid indices to local (within the current MPI
        /// process) cell indices;<br/>
        /// index: a subgrid index <em>k</em>, i.e. 0 &#8804; <em>k</em> &lt;
        /// <see cref="LocalNoOfCells"/>, or a negative number if cell
        /// <em>j</em> is not contained in the subgrid;<br/>
        /// content: a local cell index <em>j</em>, i.e. 0 &#8804;
        /// <em>j</em> &lt; <see cref="GridData.CellData.NoOfCells"/>;
        /// </summary>
        /// <remarks>
        /// The inverse of this mapping is
        /// <see cref="LocalCellIndex2SubgridIndex"/>; The array is/should be
        /// sorted in ascending order;
        /// </remarks>
        public int[] SubgridIndex2LocalCellIndex {
            get {
                if (m_SubgridIndex2LocalCellIndex == null) {
                    m_SubgridIndex2LocalCellIndex = ComputeSubgridIndex2LocalCellIndex();
                }
                return m_SubgridIndex2LocalCellIndex;
            }
        }

        /// <summary>
        /// Computes <see cref="SubgridIndex2LocalCellIndex"/>
        /// </summary>
        /// <returns></returns>
        private int[] ComputeSubgridIndex2LocalCellIndex() {
            MPICollectiveWatchDog.Watch();
            CellMask msk = m_VolumeMask;
            int[] subgridIndex2LocalCellIndex = new int[msk.NoOfItemsLocally_WithExternal];
            int jj = 0;
            foreach (Chunk c in msk.GetEnumerableWithExternal()) {
                for (int k = 0; k < c.Len; k++) {
                    subgridIndex2LocalCellIndex[jj] = c.i0 + k;
                    jj++;
                }
            }
            return subgridIndex2LocalCellIndex;
        }

        /// <summary>
        /// number of cells in the sub-grid (on current MPI process)
        /// </summary>
        public int LocalNoOfCells {
            get {
                return m_VolumeMask.NoOfItemsLocally;
            }
        }

        /// <summary>
        /// number of cells in the sub-grid (on current MPI process), including external/ghost cells
        /// </summary>
        public int LocalNoOfCells_WithExternal {
            get {
                return m_VolumeMask.NoOfItemsLocally_WithExternal;
            }
        }

        int m_GlobalNoOfCells = -1;

        /// <summary>
        /// <see cref="NoOfGhostCells"/>; computed on demand
        /// </summary>
        int m_NoOfGhostCells = -1;

        /// <summary>
        /// number of external/ghost cells
        /// </summary>
        public int NoOfGhostCells {
            get {
                if (m_GridData.iLogicalCells.NoOfExternalCells == 0)
                    // no external/ghost cells in the whole grid -> there can't be any in the subgrid
                    m_NoOfGhostCells = 0;


                if (m_NoOfGhostCells < 0) {
                    // compute number of ghost cells on demand
                    // ++++++++++++++++++++++++++++++++++++++++

                    var extMsk = m_VolumeMask.GetBitMaskWithExternal();

                    int J = m_GridData.iLogicalCells.NoOfLocalUpdatedCells;
                    int JE = extMsk.Count;
                    m_NoOfGhostCells = 0;
                    for (int j = J; j < JE; j++) {
                        if (extMsk[j])
                            m_NoOfGhostCells++;
                    }

                }

                // ret
                return m_NoOfGhostCells;
            }
        }

        /// <summary>
        /// number of cells in the sub-grid (over all MPI processes)
        /// </summary>
        public int GlobalNoOfCells {
            get {
                MPICollectiveWatchDog.Watch();

                if (m_GlobalNoOfCells < 0) {
                    m_GlobalNoOfCells = m_VolumeMask.NoOfItemsLocally.MPISum();
                }

                return m_GlobalNoOfCells;
            }
        }

        /// <summary>
        /// <see cref="VolumeMask"/>; 
        /// </summary>
        CellMask m_VolumeMask = null;

        /// <summary>
        /// all cells in the sub-grid
        /// </summary>
        public CellMask VolumeMask {
            get {
                return m_VolumeMask;
            }
        }


        EdgeMask m_BoundaryEdgesMask;

        /// <summary>
        /// all edges on the boundary of the subgrid
        /// </summary>
        public EdgeMask BoundaryEdgesMask {
            get {
                if (m_BoundaryEdgesMask == null) {

                    int E = m_GridData.iLogicalEdges.Count;
                    int[,] Edges = m_GridData.iLogicalEdges.CellIndices;

                    BitArray edges = new BitArray(E, false);
                    BitArray cells = this.VolumeMask.GetBitMaskWithExternal();

                    // loop over all Edges
                    for (int e = 0; e < E; e++) {
                        int Cel1 = Edges[e, 0];
                        int Cel2 = Edges[e, 1];

                        if (Cel2 < 0 && cells[Cel1] == true)
                            // on the domain boundary
                            edges[e] = true;
                        else {
                            if ((Cel2 >= 0) && (cells[Cel1] != cells[Cel2]))
                                edges[e] = true;
                        }
                    }

                    m_BoundaryEdgesMask = new EdgeMask(this.m_GridData, edges);
                }
                return m_BoundaryEdgesMask;
            }
        }

        EdgeMask m_InnerEdgesMask;

        /// <summary>
        /// only edges between cells in <see cref="VolumeMask"/> (see also <see cref="AllEdgesMask"/>);
        /// </summary>
        public EdgeMask InnerEdgesMask {
            get {
                if (m_InnerEdgesMask == null) {
                    int E = m_GridData.iLogicalEdges.Count;
                    int[,] Edges = m_GridData.iLogicalEdges.CellIndices;

                    BitArray edges = new BitArray(E, false);
                    BitArray cells = this.VolumeMask.GetBitMaskWithExternal();

                    // loop over all Edges
                    for (int e = 0; e < E; e++) {
                        int Cel1 = Edges[e, 0];
                        int Cel2 = Edges[e, 1];

                        if (Cel2 >= 0) {
                            if (cells[Cel1] && cells[Cel2]) {
                                edges[e] = true;
                            }
                        }
                    }

                    m_InnerEdgesMask = new EdgeMask(this.m_GridData, edges);
                }
                return m_InnerEdgesMask;
            }
        }

        
        


        EdgeMask m_AllEdgesMask;

        /// <summary>
        /// all edges with at least one cell in <see cref="VolumeMask"/>
        /// </summary>
        public EdgeMask AllEdgesMask {
            get {
                if (m_AllEdgesMask == null) {
                    m_AllEdgesMask = EdgeMask.Union(this.BoundaryEdgesMask, this.InnerEdgesMask);
                }
                return m_AllEdgesMask;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T1"></typeparam>
        /// <typeparam name="T2"></typeparam>
        /// <param name="inp_FullGridVec">
        /// input; vector associated with the "full" grid, i.e. <see cref="GridData"/>
        /// </param>
        /// <param name="outp_SubGridVec">
        /// output; vector associates with this subgrid
        /// </param>
        /// <param name="period">
        /// period, i.e.  (length/count of <paramref name="inp_FullGridVec"/>)/<see cref="LocalNoOfCells"/>
        /// </param>
        public void CompressVector<T1, T2>(T2 outp_SubGridVec, T1 inp_FullGridVec, out int period)
            where T1 : IList<double>
            where T2 : IList<double> {

            throw new NotImplementedException("todo: Christrina, wie bespr. 10feb11");
        }

        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T1"></typeparam>
        /// <param name="inp"></param>
        /// <param name="period"></param>
        /// <returns>
        /// the compressed vector, i.e. vector associated with this subgrid
        /// </returns>
        public double[] CompressVector<T1>(T1 inp, out int period)
            where T1 : IList<double> {

            throw new NotImplementedException("todo: Christrina, wie bespr. 10feb11");
        }

        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T1"></typeparam>
        /// <typeparam name="T2"></typeparam>
        /// <param name="outp_FullGridVec">
        /// output; vector associated with the "full" grid, i.e. <see cref="GridData"/>
        /// </param>
        /// <param name="inp_SubGridVec">
        /// input; vector associates with this subgrid
        /// </param>
        /// <param name="period">
        /// period, i.e.  (length/count of <paramref name="outp_FullGridVec"/>)/<see cref="LocalNoOfCells"/>
        /// </param>
        public void UncompressVector<T1, T2>(T1 outp_FullGridVec, T2 inp_SubGridVec, out int period)
            where T1 : IList<double>
            where T2 : IList<double> {

            throw new NotImplementedException("todo: Christrina, wie bespr. 10feb11");
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="RowPeriod"></param>
        /// <param name="ColPeriod"></param>
        /// <returns></returns>
        public MsrMatrix CompressOperatorMatrix(MsrMatrix inp, out int RowPeriod, out int ColPeriod) {

            throw new NotImplementedException("todo: Christrina, wie bespr. 10feb11");

        }

        /// <summary>
        /// the complement of this Subgrid
        /// </summary>
        /// <returns></returns>
        public SubGrid Complement() {
            return new SubGrid(this.VolumeMask.Complement());
        }


        double m_Volume = -1;

        /// <summary>
        /// measure of all cells within the subgrid
        /// </summary>
        public double Volume {
            get {
                if (m_Volume < 0.0) {
                    m_Volume = 0;

                    foreach (Chunk c in this.VolumeMask) {
                        int JE = c.JE;
                        for (int i = c.i0; i < JE; i++)
                            m_Volume += m_GridData.iLogicalCells.GetCellVolume(i);
                    }
                    m_Volume = m_Volume.MPISum();
                }
                return m_Volume;
            }
        }


        double m_h_min = -1;

        /// <summary>
        /// minimum cell diameter over all cells in the subgrid
        /// <see cref="GridData.CellData.h_min"/>
        /// </summary>
        public double h_minSubGrd {
            get {
                if (m_h_min < 0) {
                    var hmin = m_GridData.iGeomCells.h_min;
                    double _h_min = double.MaxValue;
                    foreach (Chunk c in this.VolumeMask) {
                        int JE = c.JE;
                        for (int i = c.i0; i < JE; i++) {
                            _h_min = Math.Min(hmin[i], _h_min);
                            ;
                        }
                    }

                    m_h_min = _h_min.MPIMin();
                }
                return m_h_min;
            }
        }

        double m_h_max = -1;

        /// <summary>
        /// maximum cell diameter over all cells in the subgrid
        /// <see cref="GridData.CellData.h_max"/>
        /// </summary>
        public double h_maxSubGrd {
            get {
                if (m_h_max < 0) {
                    var hmax = m_GridData.iGeomCells.h_max;
                    double _h_max = 0;
                    foreach (Chunk c in this.VolumeMask) {
                        int JE = c.JE;
                        for (int i = c.i0; i < JE; i++) {
                            _h_max = Math.Max(hmax[i], _h_max);
                            ;
                        }
                    }
                    m_h_max = _h_max.MPIMax();
                }
                return m_h_max;
            }
        }

        /// <summary>
        /// Tests if a matrix (associated with some operator) depends only on values in this sub-grid,
        /// or not.
        /// </summary>
        /// <param name="_Mtx">some operator matrix</param>
        /// <param name="RowMap">row-/co-domain mapping for matrix <paramref name="_Mtx"/></param>
        /// <param name="ColMap">column/domain mapping for matrix <paramref name="_Mtx"/></param>
        /// <param name="NoOfTests"></param>
        /// <returns>
        /// 0.0 if there is no dependency of the matrix <paramref name="_Mtx"/>
        /// on values outside of this subgrid.
        /// </returns>
        public double TestMatrixDependency(MsrMatrix _Mtx, UnsetteledCoordinateMapping RowMap, UnsetteledCoordinateMapping ColMap, int NoOfTests = 10) {
            Random RND = new Random();

            if (!_Mtx.RowPartitioning.Equals(RowMap))
                throw new ArgumentException();
            if (!_Mtx.ColPartition.Equals(ColMap))
                throw new ArgumentException();


            var Mtx = new ilPSP.LinSolvers.monkey.CPU.RefMatrix(_Mtx);

            int[] IR = RowMap.BasisS.Count.ForLoop(i => i);
            int[] IC = ColMap.BasisS.Count.ForLoop(i => i);

            var SgrdRow = RowMap.GetSubvectorIndices(this, true, IR);
            var SgrdCol = SgrdRow; // ColMap.GetSubvectorIndices(this, true, true, IC);

            double[] X0 = new double[ColMap.LocalLength];
            double[] Y0 = new double[RowMap.LocalLength];

            for (int i = 0; i < X0.Length; i++) {
                X0[i] = RND.NextDouble();
            }
            Mtx.SpMV(1.0, X0, 0.0, Y0);

            int i0Row = RowMap.i0;
            int i0Col = ColMap.i0;
            int iECol = ColMap.iE;
            double err = 0;
            for (int iTest = 0; iTest < NoOfTests; iTest++) {
                double[] X1 = new double[ColMap.LocalLength];
                double[] Y1 = new double[RowMap.LocalLength];

                for (int i = 0; i < X0.Length; i++) {
                    X1[i] = RND.NextDouble();
                }

                foreach (int i in SgrdCol) {
                    if (i >= i0Col && i < iECol) {
                        X1[i - i0Col] = X0[i - i0Col];
                    }
                }

                Mtx.SpMV(1.0, X1, 0.0, Y1);

                foreach (int k in SgrdRow) {
                    double Y1k = Y1[k - i0Row];
                    double Y0k = Y0[k - i0Row];

                    //Debug.Assert((Y0k - Y1k).Pow2() < 1.0e-8);

                    err += (Y0k - Y1k).Pow2();
                }
            }

            return err.MPISum();
        }
    }
}
