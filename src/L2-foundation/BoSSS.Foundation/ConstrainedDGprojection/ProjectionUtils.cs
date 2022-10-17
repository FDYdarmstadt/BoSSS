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
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP.Kraypis;

namespace BoSSS.Foundation.ConstrainedDGprojection {

    class GeometricCellForProjection {

        int CellIndex;

        int NoOfConditions;

        public GeometricCellForProjection(int cellInd) {
            CellIndex = cellInd;
            NoOfConditions = 0;
        }

        public void IncreaseNoOfConditions() {
            this.NoOfConditions++;
        }

        public int GetNoOfConditions() {
            return this.NoOfConditions;
        }

        public override bool Equals(Object obj) {

            if (obj is GeometricCellForProjection) {
                return (this.CellIndex - ((GeometricCellForProjection)obj).CellIndex) == 0;
            } else if (obj is int) {
                return (this.CellIndex - (int)obj) == 0;
            } else
                throw new ArgumentException("wrong type of object");

        }

        public override int GetHashCode() {
            return 11;
        }

    }


    class GeometricEdgeForProjection {

        public int VerticeInd1;
        public int VerticeInd2;

        //List<int> edgeList;
        int NoOfConditions;

        public GeometricEdgeForProjection(int vertInd1, int vertInd2) {
            this.VerticeInd1 = vertInd1;
            this.VerticeInd2 = vertInd2;
            //edgeList = new List<int>();
            NoOfConditions = 0;
        }

        public override bool Equals(object obj) {

            if (obj is GeometricEdgeForProjection) {

                if ((this.VerticeInd1 == ((GeometricEdgeForProjection)obj).VerticeInd1
                    && this.VerticeInd2 == ((GeometricEdgeForProjection)obj).VerticeInd2)
                    || (this.VerticeInd1 == ((GeometricEdgeForProjection)obj).VerticeInd2
                    && this.VerticeInd2 == ((GeometricEdgeForProjection)obj).VerticeInd1))
                    return true;
                else
                    return false;

            } else
                throw new ArgumentException("wrong type of object");
        }

        public override int GetHashCode() {
            return 33;
        }


        //public void AddEdge(int jEdge) {
        //    if (edgeList.Count == 4)
        //        throw new InvalidOperationException("edgeList: maximum count of 4 elements reached");
        //    edgeList.Add(jEdge);
        //}

        //public List<int> GetEdgeList() {
        //    return edgeList;
        //}

        public void IncreaseNoOfConditions() {
            this.NoOfConditions++;
        }

        public int GetNoOfConditions() {
            return this.NoOfConditions;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public int GetRefDirection(GridData m_grd, int jCell, int iFace) {

            double[] vertCoord1 = m_grd.Vertices.Coordinates.ExtractSubArrayShallow(this.VerticeInd1, -1).To1DArray();
            double[] vertCoord2 = m_grd.Vertices.Coordinates.ExtractSubArrayShallow(this.VerticeInd2, -1).To1DArray();

            int D = vertCoord1.Length;
            MultidimensionalArray vertCoord1_glb = MultidimensionalArray.Create(1, D);
            MultidimensionalArray vertCoord2_glb = MultidimensionalArray.Create(1, D);
            vertCoord1_glb.SetRow<double[]>(0, vertCoord1);
            vertCoord2_glb.SetRow<double[]>(0, vertCoord2);

            MultidimensionalArray vertCoord1_loc = MultidimensionalArray.Create(1, 1, D);
            MultidimensionalArray vertCoord2_loc = MultidimensionalArray.Create(1, 1, D);

            m_grd.TransformGlobal2Local(vertCoord1_glb, vertCoord1_loc, jCell, 1, 0);
            m_grd.TransformGlobal2Local(vertCoord2_glb, vertCoord2_loc, jCell, 1, 0);
            double[] vertCellCoord1 = vertCoord1_loc.ExtractSubArrayShallow(0, 0, -1).To1DArray();
            double[] vertCellCoord2 = vertCoord2_loc.ExtractSubArrayShallow(0, 0, -1).To1DArray();

            // identify face ref vertices
            NodeSet refV = m_grd.Edges.EdgeRefElements[0].Vertices;
            NodeSet refVvol = m_grd.Edges.EdgeRefElements[0].Vertices.GetVolumeNodeSet(m_grd, iFace, false);
            //var trafo = m_grd.Edges.Edge2CellTrafos[iFace];
            int idx1 = -1; int idx2 = -1;
            for (int i = 0; i < refV.Lengths[0]; i++) {
                double dist1 = 0.0; double dist2 = 0.0;
                for (int d = 0; d < vertCoord1.Length; d++) {
                    dist1 += (vertCellCoord1[d] - refVvol[i, d]).Pow2();
                    dist2 += (vertCellCoord2[d] - refVvol[i, d]).Pow2();
                }
                dist1 = dist1.Sqrt();
                if (dist1 < 1.0e-12)
                    idx1 = i;
                dist2 = dist2.Sqrt();
                if (dist2 < 1.0e-12)
                    idx2 = i;
            }
            Debug.Assert(idx1 != idx2);

            // identify direction of edge in ref element
            double dist = 0.0;
            for (int d = 0; d < refV.Lengths[1]; d++) {
                dist += (refV[idx1, d] - refV[idx2, d]).Pow2();
            }
            dist = dist.Sqrt();
            int dir = -1;
            for (int d = 0; d < refV.Lengths[1]; d++) {
                if (Math.Abs(Math.Abs(refV[idx1, d] - refV[idx2, d]) - dist) < 1.0e-15) {
                    dir = d;
                }
            }

            return dir;

        }

    }


    class GeometricVerticeForProjection {

        int VerticeIndex;

        int NoOfConditions;     // corresponding edge continuity constrains

        bool fixedVertice;

        double fixedValue;
        bool valueSet;

        List<int> maskedCells;
        List<int> processedBy;

        public GeometricVerticeForProjection(int vertInd) {
            VerticeIndex = vertInd;
            NoOfConditions = 0;
            fixedVertice = false;
            fixedValue = 0.0;
            valueSet = false;
            maskedCells = new List<int>();
            processedBy = new List<int>();
        }

        public int GetIndex() {
            return VerticeIndex;
        }

        public void IncreaseNoOfConditions() {
            this.NoOfConditions++;
        }

        public int GetNoOfConditions() {
            return this.NoOfConditions;
        }

        public void SetFixedVertice() {
            this.fixedVertice = true;
        }

        public bool isFixed() {
            return fixedVertice;
        }

        public void setFixedValue(double value) {
            if (valueSet) {
                fixedValue += value;
                fixedValue /= 2.0;
            } else {
                fixedValue = value;
                valueSet = true;
            }
        }

        public double getFixedValue() {
            return fixedValue;
        }

        public bool AddMaskedCell(int cellIndex) {
            if (!maskedCells.Contains(cellIndex)) {
                maskedCells.Add(cellIndex);
                return true;
            } else 
                return false;
        }

        public bool AddProcessedCell(int cellIndex) {
            if (!processedBy.Contains(cellIndex)) {
                processedBy.Add(cellIndex);
                return true;
            } else
                return false;
        }

        public bool isProcessed(int cellIndex) {
            return processedBy.Contains(cellIndex);
        }

        public List<int> GetMaskedCells() {
            return maskedCells;
        }



        public override bool Equals(Object obj) {

            if (obj is GeometricVerticeForProjection) {
                return (this.VerticeIndex - ((GeometricVerticeForProjection)obj).VerticeIndex) == 0;
            } else if (obj is int) {
                return (this.VerticeIndex - (int)obj) == 0;
            } else
                throw new ArgumentException("wrong type of object");
        }

        public override int GetHashCode() {
            return 77;
        }

    }


}
