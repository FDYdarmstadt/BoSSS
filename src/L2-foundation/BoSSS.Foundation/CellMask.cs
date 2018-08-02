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
using System.Collections;
using MPI.Wrappers;
using ilPSP.Utils;
using System.Diagnostics;
using ilPSP;
using BoSSS.Foundation.Comm;
using BoSSS.Platform;
using System.Globalization;
using System.IO;

namespace BoSSS.Foundation.Grid {


    /// <summary>
    /// masks some cells in a <see cref="GridData"/>-object
    /// </summary>
    public class CellMask : ExecutionMask {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="grddat"></param>
        /// <param name="mask">
        /// a "true" entry for all cells in grid <paramref name="grddat"/> that
        /// should be in the mask; The length of this array must not exceed
        /// <see cref="GridData.CellData.NoOfLocalUpdatedCells"/>
        /// </param>
        public CellMask(IGridData grddat, BitArray mask) :
            base(grddat, mask) {
            if (mask.Length != grddat.iLogicalCells.NoOfLocalUpdatedCells)
                throw new ArgumentException();
        }
        
        /// <summary>
        /// ctor
        /// </summary>
        public CellMask(IGridData grddat, int[] Sequence) :
            base(grddat, Sequence) 
        { }

        /// <summary>
        /// complement of this mask (all cells that are NOT in this mask);
        /// </summary>
        public CellMask Complement() {
            return base.Complement<CellMask>();
        }

        /// <summary>
        /// compiles a cell mask from a set of chunks
        /// </summary>
        /// <param name="parts">
        /// a list of chunks, which may overlap
        /// </param>
        /// <param name="grddat">
        /// the grid that this mask will be associated with;
        /// </param>
        public CellMask(IGridData grddat, params Chunk[] parts)
            : this(grddat, (IEnumerable<Chunk>)parts) {
        }

        /// <summary>
        /// compiles an quadrature execution mask from a set of chunks
        /// </summary>
        /// <param name="Parts">
        /// a list of chunks, which may overlap
        /// </param>
        /// <param name="grddat">
        /// the grid that this mask will be associated with;
        /// </param>
        protected CellMask(IGridData grddat, IEnumerable<Chunk> Parts)
            : this(grddat, FromChunkEnum(Parts))
        { }

        /// <summary>
        /// Retrieves an empty cell mask;
        /// </summary>
        /// <param name="grdDat">
        /// grid that the returned mask will be assigned to
        /// </param>
        static public CellMask GetEmptyMask(IGridData grdDat) {
            BitArray ba = new BitArray(grdDat.iLogicalCells.NoOfLocalUpdatedCells, false);
            return new CellMask(grdDat, ba);
        }

        /// <summary>
        /// Retrieves a mask containing all cells (i.e. returns the
        /// complement of <see cref="GetEmptyMask"/>)
        /// </summary>
        /// <param name="gridDat">
        /// Grid data that the returned mask will be assigned with
        /// </param>
        /// <returns>A full mask</returns>
        public static CellMask GetFullMask(IGridData gridDat) {
            return new CellMask(gridDat, new Chunk
            {
                i0 = 0,
                Len = gridDat.iLogicalCells.NoOfLocalUpdatedCells
            });
        }

        /// <summary>
        /// Selects all cells according to their cell centers, where <paramref name="SelectionFunction"/> is true
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="SelectionFunc"></param>
        /// <returns></returns>
        public static CellMask GetCellMask(Foundation.Grid.Classic.GridData gridDat, Func<double[], bool> SelectionFunc) {
            BitArray CellArray = new BitArray(gridDat.Cells.NoOfLocalUpdatedCells);
            MultidimensionalArray CellCenters = gridDat.Cells.CellCenter;
            for (int i = 0; i < gridDat.Cells.NoOfLocalUpdatedCells; i++) {
                switch (gridDat.SpatialDimension) {
                    case 1: {
                            CellArray[i] = SelectionFunc(new double[] { CellCenters[i, 0] });
                            break;
                        }
                    case 2: {
                            CellArray[i] = SelectionFunc(new double[] { CellCenters[i, 0], CellCenters[i, 1] });
                            break;
                        }
                    case 3: {
                            CellArray[i] = SelectionFunc(new double[] { CellCenters[i, 0], CellCenters[i, 1], CellCenters[i, 2] });
                            break;
                        }
                    default:
                        throw new ArgumentException();
                }
                
            }
            return new CellMask(gridDat, CellArray);
        }

        /// <summary>
        /// like the ctor.
        /// </summary>
        protected override ExecutionMask CreateInstance(IGridData grdDat, BitArray mask) {
            return new CellMask(grdDat, mask);
        }

        /// <summary>
        /// see <see cref="ExecutionMask.GetTotalNumberOfElements"/>
        /// </summary>
        protected override int GetTotalNumberOfElements(IGridData gridData) {
            return gridData.iLogicalCells.NoOfLocalUpdatedCells;
        }

        BitArray m_BitMaskWithExternal;

        /// <summary>
        /// returns a bitmask that contains also information about external/ghost cells.
        /// </summary>
        public BitArray GetBitMaskWithExternal() {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            int JE = this.GridData.iLogicalCells.Count;
            int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
            int MpiSize = this.GridData.CellPartitioning.MpiSize;

            //Debugger.Break();

            if (MpiSize > 1) {
                // more than one MPI process
                // +++++++++++++++++++++++++

                if (m_BitMaskWithExternal == null) {
                    m_BitMaskWithExternal = new BitArray(JE, false);

                    // inner cells
                    foreach (Chunk c in this) {
                        for (int i = 0; i < c.Len; i++) {
                            m_BitMaskWithExternal[i + c.i0] = true;
                        }
                    }

                    m_BitMaskWithExternal.MPIExchange(this.GridData);
                }

                return m_BitMaskWithExternal;
            } else {
                // single - processor mode
                // +++++++++++++++++++++++
                return base.GetBitMask();
            }
        }

        int m_NoOfItemsLocally_WithExternal = -1;

        /// <summary>
        /// local number of cells, including external ones
        /// </summary>
        public int NoOfItemsLocally_WithExternal {
            get {
                if (m_NoOfItemsLocally_WithExternal < 0) {
                    if (GridData.CellPartitioning.MpiSize <= 1) {
                        m_NoOfItemsLocally_WithExternal = this.NoOfItemsLocally;
                    } else {
                        this.m_NoOfItemsLocally_WithExternal = base.NoOfItemsLocally;
                        int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                        int JE = this.GridData.iLogicalCells.Count;
                        var mask = this.GetBitMaskWithExternal();

                        for (int j = J; j < JE; j++) {
                            if (mask[j])
                                this.m_NoOfItemsLocally_WithExternal++;
                        }
                    }
                }
                return m_NoOfItemsLocally_WithExternal;
            }
        }

        /// <summary>
        /// returns an enumerable structure that also contains external/ghost cells.
        /// </summary>
        public IEnumerable<Chunk> GetEnumerableWithExternal() {
            if (this.GridData.CellPartitioning.MpiSize > 1) {
                var mskExt = GetBitMaskWithExternal();
               
                var R = new List<Chunk>(this);
                int J_update = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = this.GridData.iLogicalCells.Count;
                Debug.Assert(mskExt.Count == JE);
                
                
                for (int j = J_update; j < JE; j++) {

                    if (mskExt[j]) {
                        Chunk ch;
                        ch.i0 = j;
                        ch.Len = 0;

                        while (j < JE && mskExt[j]) {
                            ch.Len++;
                            j++;
                        }

                        R.Add(ch);
                        j--;
                    }
                }

                return R;
            } else {
                return this;
            }
        }

        /// <summary>
        /// writes the center coordinates of all cells in this mask to some text file
        /// </summary>
        public override void SaveToTextFile(string fileName, bool WriteHeader = true, params ItemInfo[] infoFunc) {
            int D = GridData.SpatialDimension;
            int LI = infoFunc.Length;
            using (var file = new StreamWriter(fileName)) {
                if (WriteHeader) {

                    file.Write("Cell");
                    switch (D) {
                        case 1:
                        file.Write("\tx");
                        break;

                        case 2:
                        file.Write("\tx\ty");
                        break;

                        case 3:
                        file.Write("\tx\ty\tz");
                        break;

                        default:
                        throw new Exception();
                    }

                    for (int i = 0; i < LI; i++) {
                        file.Write("\ti(" + i + ")");
                    }
                    file.WriteLine();
                }
                double[] x = new double[D];

                foreach (int jLogCell in this.ItemEnum) { // loop over logical cells...
                    foreach (int jCell in this.GridData.GetGeometricCellIndices(jLogCell)) {
                        NodeSet localCenter = this.GridData.iGeomCells.GetRefElement(jCell).Center;
                        MultidimensionalArray globalCenters = MultidimensionalArray.Create(1, 1, D);
                        GridData.TransformLocal2Global(localCenter, jCell, 1, globalCenters, 0);

                        {
                            file.Write(jCell);
                            for (int d = 0; d < D; d++) {
                                file.Write("\t" + globalCenters[0, 0, d].ToString("e", NumberFormatInfo.InvariantInfo));
                            }

                            if (LI > 0) {
                                for (int d = 0; d < D; d++)
                                    x[d] = globalCenters[0, 0, d];
                                for (int li = 0; li < LI; li++) {
                                    double info_i = infoFunc[li](x, jLogCell, jCell);
                                    file.Write("\t" + info_i.ToString("e", NumberFormatInfo.InvariantInfo));
                                }
                            }

                            file.WriteLine();
                        }
                    }
                }
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public override bool Equals(object obj) {
            if (obj.GetType() != this.GetType())
                return false;
            return base.Equals(obj);
        }


        /// <summary>
        /// %
        /// </summary>
        public override int GetHashCode() {
            return base.GetHashCode();
        }

        /// <summary>
        /// All cells that share at least an edge with a cell in this mask.
        /// </summary>
        public CellMask AllNeighbourCells() {
            int J = this.GridData.iLogicalCells.Count;
            BitArray retMask = new BitArray(J);

            var C2E = this.GridData.iLogicalCells.Cells2Edges;
            var E2C = this.GridData.iLogicalEdges.CellIndices;

            foreach (int jCell in this.ItemEnum) {
                int[] Edges = C2E[jCell];

                foreach (int em in Edges) {
                    int iEdge = Math.Abs(em) - 1;


                    int jOtherCell, ii;
                    if (em > 0) {
                        // jCell is IN => neighbour is OUT
                        jOtherCell = E2C[iEdge, 1];
                        ii = 0;
                    } else {
                        // jCell is OUT => neighbour is IN
                        jOtherCell = E2C[iEdge, 0];
                        ii = 1;
                    }
                    Debug.Assert(E2C[iEdge, ii] == jCell);
                    
                    if(jOtherCell >= 0) // boundary edge !
                        retMask[jOtherCell] = true;
                }
            }

            return new CellMask(this.GridData, retMask);
        }

    }


}

