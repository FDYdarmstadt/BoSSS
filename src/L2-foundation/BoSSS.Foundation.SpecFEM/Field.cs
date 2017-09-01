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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;

namespace BoSSS.Foundation.SpecFEM {

    /// <summary>
    /// SpecFEM Representation of a scalar field
    /// </summary>
    public class SpecFemField {

        /// <summary>
        /// ctor;
        /// </summary>
        public SpecFemField(SpecFemBasis b) {
            m_Basis = b;
            m_Coordinates = MultidimensionalArray.Create(m_Basis.NoOfLocalNodes);
        }


        SpecFemBasis m_Basis;


        public SpecFemBasis Basis {
            get {
                return m_Basis;
            }
        }


        MultidimensionalArray m_Coordinates;

        /// <summary>
        /// Coordinates with respect to the nodal basis;<br/>
        ///  - 1st index: node index.
        /// </summary>
        public MultidimensionalArray Coordinates {
            get {
                return m_Coordinates;
            }
        }

        /// <summary>
        /// accumulate this field to a DG field
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="DGField"></param>
        /// <param name="mask">optional cell mask</param>
        public void AccToDGField(double alpha, ConventionalDGField DGField, CellMask mask = null) {
            if (!DGField.Basis.Equals(this.m_Basis.ContainingDGBasis))
                throw new ArgumentException("Basis does not match.");

            var Trafo = m_Basis.GridDat.ChefBasis.Scaling;
            var C2N = m_Basis.CellNode_To_Node;
            var MtxM2N = m_Basis.m_Modal2Nodal;
            var CellData = this.Basis.GridDat.Cells;

            int[] _K = m_Basis.NodesPerCell;
            int L = m_Basis.ContainingDGBasis.Length;
            double[][] _NodalCoordinates = _K.Select(K => new double[K]).ToArray();
            double[] ModalCoordinates = new double[L];

            if (mask == null) {
                mask = CellMask.GetFullMask(this.Basis.GridDat);
            }

            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int JE = chunk.JE;
                for (int j = j0; j < JE; j++) { // loop over cells...
                    int iKref = CellData.GetRefElementIndex(j);
                    double[] NodalCoordinates = _NodalCoordinates[iKref];
                    int K = _K[iKref];

                    // collect coordinates for cell 'j':
                    for (int k = 0; k < K; k++) {
                        int _c2n = C2N[j, k];
                        NodalCoordinates[k] = m_Coordinates[_c2n];
                    }

                    // transform
                    DGField.Coordinates.GetRow(j, ModalCoordinates);
                    MtxM2N[iKref].gemv(alpha / Trafo[j], NodalCoordinates, 1.0, ModalCoordinates);

                    // save
                    DGField.Coordinates.SetRow(j, ModalCoordinates);
                }
            }
        }


        /// <summary>
        /// projects some DG field onto this 
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="DGField"></param>
        public void ProjectDGFieldMaximum(double alpha, ConventionalDGField DGField) {

            DGField.MPIExchange();

            //var multiplicity = new int[this.m_Basis.NoOfLocalNodes];

            int J = m_Basis.GridDat.Cells.NoOfLocalUpdatedCells;
            var Trafo = m_Basis.GridDat.ChefBasis.Scaling;
            var C2N = m_Basis.CellNode_To_Node;
            var MtxN2M = m_Basis.m_Nodal2Modal;
            var CellData = this.Basis.GridDat.Cells;

            int[] _K = m_Basis.NodesPerCell;
            int L = m_Basis.ContainingDGBasis.Length;
            double[][] _NodalCoordinates = _K.Select(K => new double[K]).ToArray();
            double[] ModalCoordinates = new double[L];


            for (int j = 0; j < J; j++) { // loop over cells...
                int iKref = CellData.GetRefElementIndex(j);
                double[] NodalCoordinates = _NodalCoordinates[iKref];
                int K = _K[iKref];

                // Get DG coordinates
                Array.Clear(ModalCoordinates, 0, L);
                int Lmin = Math.Min(L, DGField.Basis.GetLength(j));
                for (int l = 0; l < Lmin; l++)
                    ModalCoordinates[l] = DGField.Coordinates[j, l];


                // transform
                DGField.Coordinates.GetRow(j, ModalCoordinates);
                MtxN2M[iKref].gemv(alpha * Trafo[j], ModalCoordinates, 0.0, NodalCoordinates);


                // collect coordinates for cell 'j':
                for (int k = 0; k < K; k++) {
                    int _c2n = C2N[j, k];
                    m_Coordinates[_c2n] = Math.Max(m_Coordinates[_c2n], NodalCoordinates[k]);
                }
            }

            using (var trx = new Transceiver(this.Basis)) {
                trx.MaxGather(m_Coordinates);
                trx.Scatter(m_Coordinates);
            }
        }


        /// <summary>
        /// projects some DG field onto this 
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="DGField"></param>
        public void ProjectDGFieldCheaply(double alpha, ConventionalDGField DGField) {

            DGField.MPIExchange();

            //var multiplicity = new int[this.m_Basis.NoOfLocalNodes];

            int J = m_Basis.GridDat.Cells.NoOfLocalUpdatedCells;
            var Trafo = m_Basis.GridDat.ChefBasis.Scaling;
            var C2N = m_Basis.CellNode_To_Node;
            var MtxN2M = m_Basis.m_Nodal2Modal;
            var CellData = this.Basis.GridDat.Cells;

            int[] _K = m_Basis.NodesPerCell;
            int L = m_Basis.ContainingDGBasis.Length;
            double[][] _NodalCoordinates = _K.Select(K => new double[K]).ToArray();
            double[] ModalCoordinates = new double[L];


            for (int j = 0; j < J; j++) { // loop over cells...
                int iKref = CellData.GetRefElementIndex(j);
                double[] NodalCoordinates = _NodalCoordinates[iKref];
                int K = _K[iKref];

                // Get DG coordinates
                Array.Clear(ModalCoordinates, 0, L);
                int Lmin = Math.Min(L, DGField.Basis.GetLength(j));
                for (int l = 0; l < Lmin; l++)
                    ModalCoordinates[l] = DGField.Coordinates[j, l];


                // transform
                DGField.Coordinates.GetRow(j, ModalCoordinates);
                MtxN2M[iKref].gemv(alpha * Trafo[j], ModalCoordinates, 0.0, NodalCoordinates);


                // collect coordinates for cell 'j':
                for (int k = 0; k < K; k++) {
                    int _c2n = C2N[j, k];
                    m_Coordinates[_c2n] += NodalCoordinates[k];
                    //multiplicity[_c2n]++;
                }
            }

            using (var trx = new Transceiver(this.Basis)) {
                trx.AccumulateGather(m_Coordinates);


                var multiplicity = this.Basis.NodeMultiplicity;
                int I = this.Basis.NoOfLocalOwnedNodes;
                for (int i = 0; i < I; i++) {
                    Debug.Assert(multiplicity[i] > 0);
                    m_Coordinates[i] *= 1.0 / multiplicity[i];

                }

                trx.Scatter(m_Coordinates);
            }
        }


        /// <summary>
        /// projects some DG field onto this
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="DGField"></param>
        /// <param name="_cm">optional restriction to computational domain</param>
        /// <remarks>
        /// This method computes an exact
        /// L2-projection of the DG-field onto the SpecFEM-space, so a global linear system, which contains all
        /// DOF's, has to be solved.
        /// In contrast, <see cref="ProjectDGFieldCheaply"/> performs an approximate projection which only involves
        /// local operations for each cell.
        /// </remarks>
        public void ProjectDGField(double alpha, ConventionalDGField DGField, CellMask _cm = null) {

            using (var trx = new Transceiver(this.Basis)) {

                CellMask cm = _cm;
                if (cm == null) {
                    cm = CellMask.GetFullMask(this.Basis.GridDat);
                }


                int J = m_Basis.GridDat.Cells.NoOfLocalUpdatedCells;
                var Trafo = m_Basis.GridDat.ChefBasis.Scaling;
                var C2N = m_Basis.CellNode_To_Node;
                var MtxM2N = m_Basis.m_Modal2Nodal;
                var CellData = this.Basis.GridDat.Cells;

                // compute RHS
                // ===========

                var b = MultidimensionalArray.Create(this.m_Basis.NoOfLocalNodes);
                {


                    int[] _K = m_Basis.NodesPerCell;
                    int L = m_Basis.ContainingDGBasis.Length;
                    double[][] _NodalCoordinates = _K.Select(K => new double[K]).ToArray(); // temporary storage for nodal coordinates per cell
                    //                                                                         1st idx: ref. elm., 2nd idx: node index
                    double[] ModalCoordinates = new double[L];


                    foreach (Chunk cnk in cm) {
                        int j0 = cnk.i0;
                        int jE = cnk.JE;
                        for (int j = j0; j < jE; j++) { // loop over cells...
                            int iKref = CellData.GetRefElementIndex(j);
                            double[] NodalCoordinates = _NodalCoordinates[iKref];
                            int K = _K[iKref];

                            if (!CellData.IsCellAffineLinear(j))
                                throw new NotSupportedException();

                            // Get DG coordinates
                            Array.Clear(ModalCoordinates, 0, L);
                            int Lmin = Math.Min(L, DGField.Basis.GetLength(j));
                            for (int l = 0; l < Lmin; l++)
                                ModalCoordinates[l] = DGField.Coordinates[j, l];

                            var tr = 1.0 / Trafo[j];

                            // transform
                            //DGField.Coordinates.GetRow(j, ModalCoordinates);
                            ModalCoordinates.ClearEntries();
                            for(int l = 0; l < Lmin; l++)
                                ModalCoordinates[l] = DGField.Coordinates[j, l];
                            MtxM2N[iKref].gemv(tr, ModalCoordinates, 0.0, NodalCoordinates, transpose: true);

                            // collect coordinates for cell 'j':
                            for (int k = 0; k < K; k++) {
                                int _c2n = C2N[j, k];
                                b[_c2n] += NodalCoordinates[k];
                            }
                        }
                    }

                }

                trx.AccumulateGather(b);

                /*

                var bcheck = new double[b.Length];
                {
                    var polys = this.Basis.NodalBasis;


                    CellQuadrature.GetQuadrature(new int[] { K },
                        this.Basis.GridDat.Context,
                        (new CellQuadratureScheme()).Compile(this.Basis.GridDat, this.Basis.ContainingDGBasis.Degree*2),
                        delegate(MultidimensionalArray NodesUntransformed) { // Del_CreateNodeSetFamily
                            var NSC = this.Basis.GridDat.Context.NSC;
                            return new NodeSetController.NodeSetContainer[] { NSC.CreateContainer(NodesUntransformed) };
                        },
                        delegate(int i0, int Length, int _NoOfNodes, MultidimensionalArray EvalResult) {
                            var PolyAtNode = MultidimensionalArray.Create(K, _NoOfNodes);
                            for (int k = 0; k < K; k++) {
                                polys[k].Evaluate(PolyAtNode.ExtractSubArrayShallow(k, -1), this.Basis.GridDat.Context.NSC.Current_NodeSetFamily[0].NodeSet);
                            }

                            var DGFatNodes = MultidimensionalArray.Create(Length, _NoOfNodes);
                            DGField.Evaluate(i0, Length, 0, DGFatNodes);

                            //for(int i = 0; i < Length; i++) {
                            //    for (int n = 0; n < _NoOfNodes; n++) {
                            //        for (int k = 0; k < K; k++) {
                            //            for (int l = 0; l < K; l++) {
                            //                EvalResult[i, n, k, l] = PolyAtNode[k, n]*PolyAtNode[l, n]; 
                            //            }
                            //        }
                            //    }
                            //}

                            EvalResult.Multiply(1.0, PolyAtNode, DGFatNodes, 0.0, "jnk", "kn", "jn");

                            //double errSum = 0;
                            //for (int i = 0; i < Length; i++) {
                            //    for (int n = 0; n < _NoOfNodes; n++) {
                            //        for (int k = 0; k < K; k++) {
                            //            for (int l = 0; l < K; l++) {
                            //                double soll = PolyAtNode[k, n]*PolyAtNode[l, n];
                            //                errSum += Math.Abs(soll - EvalResult[i, n, k, l]);
                            //            }
                            //        }
                            //    }
                            //}
                            //Console.WriteLine("errsum = " + errSum);
                        },
                        delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                            for (int i = 0; i < Length; i++) {
                                int jCell = i + i0;

                                for (int k = 0; k < K; k++) {
                                    bcheck[C2N[jCell, k]] += ResultsOfIntegration[i, k];
                                }
                            
                                //CellMass[jCell] = new FullMatrix(K, K);
                                //CellMass[jCell].Initialize(ResultsOfIntegration.ExtractSubArrayShallow(i, -1, -1));
                            }
                        }).Execute();


                    double fuck = GenericBlas.L2Dist(b, bcheck);
                    Console.WriteLine("Distance error = " + fuck);

                }


                */

                if (_cm == null) {
                    // full domain projection branch
                    // +++++++++++++++++++++++++++++


                    var x = new double[this.Basis.NoOfLocalOwnedNodes];
                    var solStat = m_Basis.MassSolver.Solve(x, b.ExtractSubArrayShallow(new int[] { 0 }, new int[] { this.Basis.NoOfLocalOwnedNodes - 1 }).To1DArray());

                    {
                        if (solStat.Converged == false)
                            throw new ArithmeticException("DG -> SpecFEM Projection failed because the Mass matrix solver did not converge.");


                        double[] chk = b.ExtractSubArrayShallow(new int[] { 0 }, new int[] { this.Basis.NoOfLocalOwnedNodes - 1 }).To1DArray();
                        this.Basis.MassMatrix.SpMVpara(-1.0, x, 1.0, chk);
                        double chk_nomr = chk.L2Norm();

                        if (chk_nomr >= 1.0e-8) {
                            throw new ArithmeticException(string.Format("DG -> SpecFEM Projection failed: solver converged, but with high residual {0}.", chk_nomr.ToStringDot()));
                        }

                    }

                    //m_Basis.MassMatrix.SpMV(1.0, b, 0.0, x);
                    m_Coordinates.ExtractSubArrayShallow(new int[] { 0 }, new int[] { this.Basis.NoOfLocalOwnedNodes - 1 }).AccVector(alpha, x);
                    //m_Coordinates.AccV(alpha, b);

                } else {
                    // restricted domain projection branch
                    // +++++++++++++++++++++++++++++++++++

                    List<int> OccupiedRows_Global = new List<int>();
                    //List<int> OccupiedRows_Local = new List<int>();

                    var MM = Basis.ComputeMassMatrix(cm);
                    int i0 = MM.RowPartitioning.i0, iE = MM.RowPartitioning.iE;
                    for (int i = i0; i < iE; i++) {
                        if (MM.GetNoOfNonZerosPerRow(i) > 0) {
                            OccupiedRows_Global.Add(i);
                            //OccupiedRows_Local.Add(i - i0);
                        }
                    }

                    var CompressedPart = new Partitioning(OccupiedRows_Global.Count);
                    var CompressedMM = new MsrMatrix(CompressedPart);

                    MM.WriteSubMatrixTo(CompressedMM, OccupiedRows_Global, default(int[]), OccupiedRows_Global, default(int[]));

                    var b_sub = new double[OccupiedRows_Global.Count];
                    //try {
                        b_sub.AccV(1.0, b.To1DArray(), default(int[]), OccupiedRows_Global, b_index_shift: -i0);
                    //} catch(Exception e) {
                    //    Debugger.Launch();
                    //}
                    //csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                    var x_sub = new double[b_sub.Length];

                    var solver = new ilPSP.LinSolvers.monkey.CG();
                    solver.MatrixType = ilPSP.LinSolvers.monkey.MatrixType.CCBCSR;
                    solver.DevType = ilPSP.LinSolvers.monkey.DeviceType.CPU;
                    solver.ConvergenceType = ConvergenceTypes.Absolute;
                    solver.Tolerance = 1.0e-12;
                    solver.DefineMatrix(CompressedMM);

                    var solStat = solver.Solve(x_sub, b_sub.CloneAs());
                    {
                        if (solStat.Converged == false)
                            throw new ArithmeticException("DG -> SpecFEM Projection failed because the Mass matrix solver did not converge.");

                        var chk = b_sub;
                        CompressedMM.SpMVpara(-1.0, x_sub, 1.0, chk);
                        double chk_nomr = chk.L2Norm();

                        if (chk_nomr >= 1.0e-8) {
                            throw new ArithmeticException(string.Format("DG -> SpecFEM Projection failed: solver converged, but with high residual {0}.", chk_nomr.ToStringDot()));
                        }
                    }

                    double[] x = new double[this.Basis.NoOfLocalOwnedNodes];
                    x.AccV(1.0, x_sub, OccupiedRows_Global, default(int[]), acc_index_shift: -i0);
                    m_Coordinates.ExtractSubArrayShallow(new int[] { 0 }, new int[] { this.Basis.NoOfLocalOwnedNodes - 1 }).AccVector(alpha, x);
                }

                trx.Scatter(this.m_Coordinates);
            }


        }
    }
}
