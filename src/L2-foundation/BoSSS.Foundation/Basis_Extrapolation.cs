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
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Foundation {

    public partial class Basis {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="CellPairs">
        /// 1st index: list of cells <br/>
        /// 2nd index: in {0, 1}
        /// </param>
        /// <param name="M">
        /// 1st index: corresponds with 1st index of <paramref name="CellPairs"/><br/>
        /// 2nd index: matrix row index <br/>
        /// 3rd index: matrix column index
        /// </param>
        /// <param name="Minv">the inverse of <paramref name="M"/></param>
        /// <remarks>
        /// Let \f$ K_j \f$ and \f$ K_i \f$ be two different cells with a linear-affine 
        /// transformation to the reference element.
        /// Here, \f$ j \f$=<paramref name="CellPairs"/>[a,0] and \f$ i \f$=<paramref name="CellPairs"/>[a,1]. 
        /// The DG-basis in these cells can uniquely be represented as
        /// \f[ 
        /// \phi_{j n} (\vec{x}) = p_n (\vec{x}) \vec{1}_{K_j} (\vec{x})
        /// \textrm{ and }
        /// \phi_{i m} (\vec{x}) = q_m (\vec{x}) \vec{1}_{K_i} (\vec{x})
        /// \f]
        /// where \f$ \vec{1}_X \f$ denotes the characteristic function for set \f$ X \f$
        /// and \f$ p_n\f$  and \f$ p_m\f$  are polynomials.
        /// Then, for the output \f$ M \f$ =<paramref name="M"/>[a,-,-] fulfills
        /// \f[ 
        /// \phi_{j n} + \sum_{m} M_{m n} \phi_{i m}
        /// =
        /// p_n \vec{1}_{K_j \cup K_i}
        /// \f]
        /// </remarks>
        public void GetExtrapolationMatrices(int[,] CellPairs, MultidimensionalArray M, MultidimensionalArray Minv = null) {
            var m_Context = this.GridDat;
            int N = this.Length;
            int Esub = CellPairs.GetLength(0);
            int JE = this.GridDat.iLogicalCells.Count;
            int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

            if (CellPairs.GetLength(1) != 2)
                throw new ArgumentOutOfRangeException("second dimension is expected to be 2!");
            if (M.Dimension != 3)
                throw new ArgumentException();
            if (M.GetLength(0) != Esub)
                throw new ArgumentException();
            if (M.GetLength(1) != N || M.GetLength(2) != N)
                throw new ArgumentException();
            if (Minv != null) {
                if (Minv.GetLength(0) != Esub)
                    throw new ArgumentException();
                if (Minv.GetLength(1) != N || Minv.GetLength(2) != N)
                    throw new ArgumentException();
            }

            MultidimensionalArray NodesGlobal = new MultidimensionalArray(3);
            MultidimensionalArray Minv_tmp = MultidimensionalArray.Create(N, N);
            MultidimensionalArray M_tmp = MultidimensionalArray.Create(N, N);


            for (int esub = 0; esub < Esub; esub++) { // loop over the cell pairs...

                int jCell0 = CellPairs[esub, 0];
                int jCell1 = CellPairs[esub, 1];
                if (jCell0 < 0 || jCell0 >= JE)
                    throw new ArgumentOutOfRangeException("Cell index out of range.");
                if (jCell1 < 0 || jCell1 >= JE)
                    throw new ArgumentOutOfRangeException("Cell index out of range.");

                bool swap;
                if (jCell0 >= J) {
                    //if(true) {
                    swap = true;
                    int a = jCell0;
                    jCell0 = jCell1;
                    jCell1 = a;
                } else {
                    swap = false;
                }


                if (!m_Context.iGeomCells.IsCellAffineLinear(jCell0))
                    throw new NotSupportedException("Currently not supported for curved cells.");
                if (!m_Context.iGeomCells.IsCellAffineLinear(jCell1))
                    throw new NotSupportedException("Currently not supported for curved cells.");

                Debug.Assert(jCell0 < J);

                var cellMask = new CellMask(m_Context, new[] { new Chunk() { i0 = jCell0, Len = 1 } }, MaskType.Geometrical);

                // we project the basis function from 'jCell1' onto 'jCell0'

                CellQuadrature.GetQuadrature(new int[2] { N, N }, m_Context,
                    (new CellQuadratureScheme(true, cellMask)).Compile(m_Context, this.Degree * 2), // integrate over target cell
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                        NodeSet nodes_Cell0 = QR.Nodes;
                        Debug.Assert(Length == 1);

                        NodesGlobal.Allocate(1, nodes_Cell0.GetLength(0), nodes_Cell0.GetLength(1));
                        m_Context.TransformLocal2Global(nodes_Cell0, jCell0, 1, NodesGlobal, 0);
                        var nodes_Cell1 = new NodeSet(GridDat.iGeomCells.GetRefElement(jCell1), nodes_Cell0.GetLength(0), nodes_Cell0.GetLength(1));
                        m_Context.TransformGlobal2Local(NodesGlobal.ExtractSubArrayShallow(0, -1, -1), nodes_Cell1, jCell1, null);
                        nodes_Cell1.LockForever();


                        var phi_0 = this.CellEval(nodes_Cell0, jCell0, 1).ExtractSubArrayShallow(0, -1, -1);
                        var phi_1 = this.CellEval(nodes_Cell1, jCell1, 1).ExtractSubArrayShallow(0, -1, -1);

                        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1);

                        EvalResult.Multiply(1.0, phi_1, phi_0, 0.0, "kmn", "kn", "km");
                    },
                    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                     Debug.Assert(Length == 1);

                                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1);
                                                     Minv_tmp.Clear();
                                                     Minv_tmp.Acc(1.0, res);
                                                 }).Execute();

                // compute the inverse
                Minv_tmp.InvertTo(M_tmp);

                // store
                if (!swap) {
                    M.ExtractSubArrayShallow(esub, -1, -1).AccMatrix(1.0, M_tmp);
                    if (Minv != null) {
                        Minv.ExtractSubArrayShallow(esub, -1, -1).AccMatrix(1.0, Minv_tmp);
                    }
                } else {
                    M.ExtractSubArrayShallow(esub, -1, -1).AccMatrix(1.0, Minv_tmp);
                    if (Minv != null) {
                        Minv.ExtractSubArrayShallow(esub, -1, -1).AccMatrix(1.0, M_tmp);
                    }
                }
            }
        }

        /// <summary>
        /// Change-of-basis, in cell 0 
        /// </summary>
        /// <param name="jCell"></param>
        /// <param name="pl"></param>
        public MultidimensionalArray GetChangeofBasisMatrix(int jCell, PolynomialList pl) {
            var m_Context = this.GridDat;
            int N = this.Length;
            int M = pl.Count;
            int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

            if (jCell < 0 || jCell >= J)
                throw new ArgumentOutOfRangeException("cell index out of range");
            MultidimensionalArray Mtx = MultidimensionalArray.Create(N, M);
            
            
            var cellMask = new CellMask(m_Context, new[] { new Chunk() { i0 = jCell, Len = 1 } }, MaskType.Geometrical);

            // we project the basis function from 'jCell1' onto 'jCell0'

            CellQuadrature.GetQuadrature(new int[2] { N, M }, m_Context,
                (new CellQuadratureScheme(true, cellMask)).Compile(m_Context, this.Degree * 2), // integrate over target cell
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                    NodeSet nodes_Cell0 = QR.Nodes;
                    Debug.Assert(Length == 1);

                    //NodesGlobal.Allocate(1, nodes_Cell0.GetLength(0), nodes_Cell0.GetLength(1));
                    //m_Context.TransformLocal2Global(nodes_Cell0, jCell0, 1, NodesGlobal, 0);
                    //var nodes_Cell1 = new NodeSet(GridDat.iGeomCells.GetRefElement(jCell1), nodes_Cell0.GetLength(0), nodes_Cell0.GetLength(1));
                    //m_Context.TransformGlobal2Local(NodesGlobal.ExtractSubArrayShallow(0, -1, -1), nodes_Cell1, jCell1, null);
                    //nodes_Cell1.LockForever();

                    var phi_0 = this.CellEval(nodes_Cell0, jCell, 1).ExtractSubArrayShallow(0, -1, -1);
                    MultidimensionalArray R = MultidimensionalArray.Create(QR.NoOfNodes, pl.Count);
                    pl.Evaluate(nodes_Cell0, R);

                    var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1);
                    EvalResult.Multiply(1.0, R, phi_0, 0.0, "knm", "km", "kn");
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    Debug.Assert(Length == 1);

                    var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1);
                    Mtx.Clear();
                    Mtx.Acc(1.0, res);
                }).Execute();

            return Mtx;
        }

    }
}
