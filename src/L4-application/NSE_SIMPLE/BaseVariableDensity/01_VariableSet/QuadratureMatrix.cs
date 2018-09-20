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
using BoSSS.Foundation;
using BoSSS.Platform;
using BoSSS.Foundation.Quadrature;
using ilPSP.LinSolvers;
using BoSSS.Foundation.Grid;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace NSE_SIMPLE {

    /// <summary>
    /// Matrix representation of a DG Field
    /// or of a function of this Field.
    /// </summary>
    public class QuadratureMatrix : CellQuadrature {

        Basis m_Basis;

        SinglePhaseField[] m_Fields;

        BlockDiagonalMatrix m_Matrix;

        public BlockDiagonalMatrix Matrix {
            get {
                return m_Matrix;
            }
        }

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="Basis"></param>        
        /// <param name="GridDat"></param>
        /// <param name="Fields"></param>
        public QuadratureMatrix(Basis Basis, IGridData GridDat, params SinglePhaseField[] Fields)
            : base(new int[] { Basis.MaximalLength * Basis.MaximalLength },
                GridDat,
                (new CellQuadratureScheme()).Compile(GridDat, Fields[0].Basis.Degree + 2 * Basis.Degree)) //
        {
            m_Basis = Basis;
            m_Fields = Fields;

            UnsetteledCoordinateMapping map = new UnsetteledCoordinateMapping(Basis);
            m_Matrix = new BlockDiagonalMatrix(map.LocalLength, Basis.MaximalLength);

            FieldVals = new MultidimensionalArray[Fields.Length];
        }

        MultidimensionalArray BasisValues;

        protected override void QuadNodesChanged(NodeSet newNodes) {
            BasisValues = m_Basis.Evaluate(newNodes);
        }

        MultidimensionalArray[] FieldVals;

        protected override void AllocateBuffers(int NoOfItems, NodeSet ruleNodes) {
            base.AllocateBuffers(NoOfItems, ruleNodes);
            for (int i = 0; i < FieldVals.Length; i++) {
                FieldVals[i] = MultidimensionalArray.Create(base.gridData.iLogicalCells.NoOfLocalUpdatedCells, ruleNodes.GetLength(0));
            }
        }

        
        protected override void Evaluate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {

            int N = m_Basis.MaximalLength;
            int NoOfNodes = QR.Nodes.NoOfNodes;

            if (((GridData)gridData).Cells.ContainsNonlinearCell())
                throw new NotImplementedException();

            MultidimensionalArray scaling = base.gridData.iGeomCells.JacobiDet;

            for (int i = 0; i < m_Fields.Length; i++) {
                m_Fields[i].Evaluate(i0, Length, QR.Nodes, FieldVals[i], 0, 0.0);
            }

            for (int j = 0; j < Length; j++) {
                int jCell = j + i0;
                double sc = 1.0 / scaling[jCell];

                for (int node = 0; node < NoOfNodes; node++) {
                    for (int BlkRowLoc = 0; BlkRowLoc < N; BlkRowLoc++) {
                        for (int BlkColLoc = 0; BlkColLoc < N; BlkColLoc++) {
                            double[] _FieldVals = new double[FieldVals.Length];
                            for (int i = 0; i < FieldVals.Length; i++) {
                                _FieldVals[i] = FieldVals[i][j, node];
                            }
                            double FieldFuncValue = FieldFunc(_FieldVals);
                            EvalResult[j, node, BlkRowLoc * N + BlkColLoc] = BasisValues[node, BlkRowLoc] * BasisValues[node, BlkColLoc] * sc * FieldFuncValue;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Identity function - can be overwritten in derived classes.
        /// </summary>
        /// <param name="FieldVal"></param>
        /// <returns></returns>
        protected virtual double FieldFunc(params double[] FieldVal) {
            return FieldVal[0];
        }

        protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {

            int N = m_Basis.MaximalLength;
            int i0_matrix = (int)m_Matrix.RowPartitioning.i0;

            for (int j = 0; j < Length; j++) {
                int jCell = j + i0;

                for (int BlkRowLoc = 0; BlkRowLoc < N; BlkRowLoc++) {
                    int MtxRowGlobal = i0_matrix + jCell * N + BlkRowLoc;
                    for (int BlkColLoc = 0; BlkColLoc < N; BlkColLoc++) {
                        int MtxColGlobal = i0_matrix + jCell * N + BlkColLoc;
                        m_Matrix[MtxRowGlobal, MtxColGlobal] = ResultsOfIntegration[j, BlkRowLoc * N + BlkColLoc];
                    }
                }
            }
        }

        public void Update() {
            this.Execute();
        }
    }
}