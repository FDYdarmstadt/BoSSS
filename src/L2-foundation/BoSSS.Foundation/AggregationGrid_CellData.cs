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
using System.Threading.Tasks;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Aggregation {
    partial class AggregationGrid {

        public IGeometricalCellsData iGeomCells {
            get {
                return m_GeomCellData;
            }
        }

        GeomCellData m_GeomCellData;

        /// <summary>
        /// Just a wrapper/proxy around the geometrical cell data (<see cref="IGridData.iGeomCells"/>) 
        /// of the parent grid (<see cref="AggregationGrid.ParentGrid"/>).
        /// </summary>
        class GeomCellData : IGeometricalCellsData {
            internal AggregationGrid m_Owner;

            public int[][] CellVertices {
                get {
                    return m_Owner.ParentGrid.iGeomCells.CellVertices;
                }
            }

            public MultidimensionalArray h_max {
                get {
                    return m_Owner.ParentGrid.iGeomCells.h_max;
                }
            }

            public MultidimensionalArray h_min {
                get {
                    return m_Owner.ParentGrid.iGeomCells.h_min;
                }
            }

            public CellInfo[] InfoFlags {
                get {
                    return m_Owner.ParentGrid.iGeomCells.InfoFlags;
                }
            }

            public MultidimensionalArray InverseTransformation {
                get {
                    return m_Owner.ParentGrid.iGeomCells.InverseTransformation;
                }
            }

            public MultidimensionalArray JacobiDet {
                get {
                    return m_Owner.ParentGrid.iGeomCells.JacobiDet;
                }
            }

            public int Count {
                get {
                    return m_Owner.ParentGrid.iGeomCells.Count;
                }
            }

            public RefElement[] RefElements {
                get {
                    return m_Owner.ParentGrid.iGeomCells.RefElements;
                }
            }

            public MultidimensionalArray Transformation {
                get {
                    return m_Owner.ParentGrid.iGeomCells.Transformation;
                }
            }

            public int[] GeomCell2LogicalCell {
                get;
                internal set;
            }

            public void GetCellBoundingBox(int j, BoundingBox bb) {
                m_Owner.ParentGrid.iGeomCells.GetCellBoundingBox(j, bb);
            }


            CellMask[] m_Cells4Refelement;

            public CellMask GetCells4Refelement(RefElement Kref) {
                int iKref = Array.IndexOf(this.RefElements, Kref);

                if (m_Cells4Refelement == null) {
                    m_Cells4Refelement = new CellMask[this.RefElements.Length];
                }

                if(m_Cells4Refelement[iKref] == null) { 
                    var OrgMsk = m_Owner.ParentGrid.iGeomCells.GetCells4Refelement(Kref);
                    Debug.Assert(OrgMsk.MaskType == MaskType.Geometrical);
                    Debug.Assert(object.ReferenceEquals(OrgMsk.GridData, m_Owner.ParentGrid));

                    m_Cells4Refelement[iKref] = new CellMask(m_Owner, OrgMsk, MaskType.Geometrical);
                }
                return m_Cells4Refelement[iKref];
            }

            public CellType GetCellType(int jCell) {
                return m_Owner.ParentGrid.iGeomCells.GetCellType(jCell);
            }

            public double GetCellVolume(int j) {
                return m_Owner.ParentGrid.iGeomCells.GetCellVolume(j);
            }

            public int GetNoOfSimilarConsecutiveCells(CellInfo mask, int j0, int Lmax) {
                throw new NotImplementedException();
            }

            public RefElement GetRefElement(int j) {
                return this.RefElements[GetRefElementIndex(j)];
            }

            public int GetRefElementIndex(int jCell) {
                return m_Owner.ParentGrid.iGeomCells.GetRefElementIndex(jCell);
            }

            public bool IsCellAffineLinear(int j) {
                return m_Owner.ParentGrid.iGeomCells.IsCellAffineLinear(j);
            }
        }

        LogicalCellData m_LogicalCellData;

        public ILogicalCellData iLogicalCells {
            get {
                return m_LogicalCellData;
            }
        }

        public class LogicalCellData : ILogicalCellData {
            internal AggregationGrid m_Owner;

            public int[][] AggregateCellToParts {
                get;
                internal set;
            }

            public int[][] CellNeighbours {
                get;
                internal set;
            }

            public int[][] Cells2Edges {
                get;
                internal set;
            }

            public int Count {
                get {
                    return NoOfLocalUpdatedCells + NoOfExternalCells;
                }
            }

            public int NoOfExternalCells {
                get {
                    return m_Owner.m_Parallel.GlobalIndicesExternalCells.Length;
                }
            }

            public int NoOfLocalUpdatedCells {
                get {
                    Debug.Assert(m_Owner.CellPartitioning.LocalLength == CellNeighbours.Length);
                    return m_Owner.CellPartitioning.LocalLength;
                }
            }

            public void GetCellBoundingBox(int j, BoundingBox bb) {
                BoundingBox bbTemp = new BoundingBox(bb.D);
                bb.Clear();
                foreach(int jPart in this.AggregateCellToParts[j]) {
                    m_Owner.m_GeomCellData.GetCellBoundingBox(jPart, bbTemp);
                    bb.AddBB(bbTemp);
                }
            }

            public CellMask GetCells4Refelement(RefElement Kref) {
                throw new NotImplementedException();
            }

            public double GetCellVolume(int j) {
                double Sum = 0;
                foreach(int jPart in AggregateCellToParts[j]) {
                    Sum += m_Owner.m_GeomCellData.GetCellVolume(j);
                }
                return Sum;
            }

            public long GetGlobalID(int j) {
                throw new NotImplementedException();
            }

            public bool IsCellAffineLinear(int j) {
                bool ret = true;
                foreach(int jPart in AggregateCellToParts[j]) {
                    ret &= m_Owner.m_GeomCellData.IsCellAffineLinear(j);
                }
                return ret;
            }
        }

    }
}
