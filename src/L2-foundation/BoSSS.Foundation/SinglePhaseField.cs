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
using System.Runtime.Serialization;
using BoSSS.Foundation.IO;
using BoSSS.Platform;
using ilPSP.Tracing;
using ilPSP.Utils;
using ilPSP;

namespace BoSSS.Foundation {

    /// <summary>
    /// 'normal' DG field
    /// </summary>
    public class SinglePhaseField : ConventionalDGField {

        /// <summary>
        /// an implementation of <see cref="FieldFactory{T}"/> that creates
        /// <see cref="SinglePhaseField"/>-DG-fields.
        /// </summary>
        /// <param name="__Basis">The basis that is used for this field</param>
        /// <param name="__Identification">
        /// identification string for this field; This can be null or empty,
        /// however, if IO should be performed for this object, the
        /// identification must be unique within the given context.
        /// </param>
        /// <returns>a <see cref="SinglePhaseField"/>-instance</returns>
        public static SinglePhaseField Factory(Basis __Basis, String __Identification) {
            return new SinglePhaseField(__Basis, __Identification);
        }

        /// <summary>
        /// constructs a new field with an empty identification string 
        /// (see <see cref="DGField.Identification"/>); Because of this, this
        /// field can't be used for IO;
        /// </summary>
        public SinglePhaseField(Basis __Basis) :
            this(__Basis, null) {
        }

        /// <summary>
        /// storage for DG coordinates: <see cref="DGField.Coordinates"/>
        /// </summary>
        [NonSerialized]
        MultidimensionalArray m_Mda_Coordinates = new MultidimensionalArray(2);

        /// <summary>
        /// DG Coordinates
        /// </summary>
        public override IMatrix Coordinates {
            get {
                return m_Mda_Coordinates;
            }
        }

        /// <summary>
        /// see <see cref="DGField.Evaluate(int,int,NodeSet,MultidimensionalArray,int,double)"/>
        /// </summary>
        override public void Evaluate(int j0, int Len, NodeSet NS, MultidimensionalArray result, int ResultIndexOffset, double ResultPreScale) {
            int D = GridDat.SpatialDimension; // spatial dimension
            int N = m_Basis.Length;      // number of coordinates per cell
            int M = NS.NoOfNodes;        // number of nodes

            if(result.Dimension != 2)
                throw new ArgumentOutOfRangeException("result", "dimension of result array must be 2");
            if(Len > result.GetLength(0) + ResultIndexOffset)
                throw new ArgumentOutOfRangeException("mismatch between Len and 0-th length of result");
            if(result.GetLength(1) != M)
                throw new ArgumentOutOfRangeException();

            var coöSys = NS.GetNodeCoordinateSystem(this.GridDat);
            var resultAcc = result.ExtractSubArrayShallow(new int[] { ResultIndexOffset, 0 }, new int[] { ResultIndexOffset + Len - 1, M - 1 });

            if(coöSys != NodeCoordinateSystem.CellCoord)
                throw new ArgumentOutOfRangeException();

            var coord = m_Mda_Coordinates.ExtractSubArrayShallow(new int[] { j0, 0 }, new int[] { j0 + Len - 1, N - 1 });
            DGField.EvaluateInternal(j0, Len, NS, this.m_Basis, coord, resultAcc, ResultPreScale);


        }

        /// <summary>
        /// Evaluates the field along edges.
        /// </summary>
        /// <param name="ResultIndexOffset">
        /// An offset for the first index of <paramref name="ValueIN"/> resp. <paramref name="ValueOT"/>,
        /// i.e. the first result will be written to
        /// <paramref name="ValueIN"/>[<paramref name="ResultIndexOffset"/>,*].
        /// </param>
        /// <param name="e0">Index of the first edge to evaluate.</param>
        /// <param name="Len">Number of edges to evaluate</param>
        /// <param name="NS">
        /// nodes to evaluate at
        /// </param>
        /// <param name="ValueIN">
        /// If not null, contains the following output: 
        /// On exit, the value of the DG field at the given nodes are
        /// <em>accumulated</em> (!) there. Before the values are added, 
        /// the original content is scaled by <paramref name="ResultPreScale"/>.<br/>
        /// The array is 2-dimensional:
        /// <list type="bullet">
        ///   <item>
        ///   1st index: edge index <i>j</i> - <paramref name="e0"/>;
        ///   </item>
        ///   <item>
        ///   2nd index: node index <i>k</i>, corresponds with 1st index of
        ///   the node set <paramref name="NS"/>;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="ValueOUT">
        /// Same as <paramref name="ValueIN"/>.
        /// </param>
        /// <param name="MeanValueIN">
        /// If not null, contains the following output: 
        /// On exit, the mean values of the DG field at the given edges are
        /// <em>accumulated</em> (!) there. Before the values are added, 
        /// the original content is scaled by <paramref name="ResultPreScale"/>.<br/>
        /// The array is 2-dimensional:
        /// <list type="bullet">
        ///   <item>
        ///   1st index: edge index <i>j</i> - <paramref name="e0"/>;
        ///   </item>
        ///   <item>
        ///   2nd index: node index <i>k</i>, corresponds with 1st index of
        ///   the node set <paramref name="NS"/>;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="MeanValueOT">
        /// Same as <paramref name="MeanValueIN"/>.
        /// </param>
        /// <param name="GradientIN">
        /// If not null, contains the following output: 
        /// On exit, the value of the gradient of DG field at the given nodes
        /// are <em>accumulated</em> (!) there. Before the values are added, 
        /// the original content is scaled by <paramref name="ResultPreScale"/>.<br/>
        /// The array is 3-dimensional:
        /// <list type="bullet">
        ///   <item>
        ///   1st index: edge index <i>j</i> - <paramref name="e0"/>;
        ///   </item>
        ///   <item>
        ///   2nd index: node index <i>k</i>, corresponds with 1st index of
        ///   the node set <paramref name="NS"/>;
        ///   </item>
        ///   <item>
        ///   2rd index: spatial dimension;
        ///   </item>
        /// </list>
        /// </param>
        /// <param name="GradientOT">
        /// Same as <paramref name="GradientIN"/>.
        /// </param>
        /// <param name="ResultPreScale">
        /// Scaling that is applied to <paramref name="ValueIN"/> and <paramref name="ValueOUT"/> before
        /// the field evaluation is added
        /// </param>
        public override void EvaluateEdge(int e0, int Len, NodeSet NS,
            MultidimensionalArray ValueIN, MultidimensionalArray ValueOT,
            MultidimensionalArray MeanValueIN = null, MultidimensionalArray MeanValueOT = null,
            MultidimensionalArray GradientIN = null, MultidimensionalArray GradientOT = null,
            int ResultIndexOffset = 0, double ResultPreScale = 0.0) //
        {


            int D = GridDat.SpatialDimension; // spatial dimension
            int N = m_Basis.Length;      // number of coordinates per cell
            int K = NS.NoOfNodes;        // number of nodes

            if((ValueIN != null) != (ValueOT != null))
                throw new ArgumentException();
            if((MeanValueIN != null) != (MeanValueOT != null))
                throw new ArgumentException();
            if((GradientIN != null) != (GradientOT != null))
                throw new ArgumentException();

            if(ValueIN != null) {
                if(ValueIN.Dimension != 2)
                    throw new ArgumentException("result", "dimension of result array must be 2");
                if(Len > ValueIN.GetLength(0) + ResultIndexOffset)
                    throw new ArgumentException("mismatch between Len and 0-th length of result");
                if(ValueIN.GetLength(1) != K)
                    throw new ArgumentException();
                if(ValueOT.Dimension != 2)
                    throw new ArgumentException("result", "dimension of result array must be 2");
                if(Len > ValueOT.GetLength(0) + ResultIndexOffset)
                    throw new ArgumentException("mismatch between Len and 0-th length of result");
                if(ValueOT.GetLength(1) != K)
                    throw new ArgumentException();

                if(!(ResultIndexOffset == 0 && Len == ValueIN.GetLength(0)))
                    ValueIN = ValueIN.ExtractSubArrayShallow(new int[] { ResultIndexOffset, 0 }, new int[] { ResultIndexOffset + Len - 1, K - 1 });

                if(!(ResultIndexOffset == 0 && Len == ValueOT.GetLength(0)))
                    ValueOT = ValueOT.ExtractSubArrayShallow(new int[] { ResultIndexOffset, 0 }, new int[] { ResultIndexOffset + Len - 1, K - 1 });
            }

            if(MeanValueIN != null) {
                if(MeanValueIN.Dimension != 1)
                    throw new ArgumentException("Dimension of mean-value result array must be 1.");
                if(Len > MeanValueIN.GetLength(0) + ResultIndexOffset)
                    throw new ArgumentException("mismatch between Len and 0-th length of result");
                if(MeanValueOT.Dimension != 1)
                    throw new ArgumentException("Dimension of mean-value result array must be 1.");
                if(Len > MeanValueOT.GetLength(0) + ResultIndexOffset)
                    throw new ArgumentException("mismatch between Len and 0-th length of result");

                if(!(ResultIndexOffset == 0 && Len == MeanValueIN.GetLength(0)))
                    MeanValueIN = MeanValueIN.ExtractSubArrayShallow(new int[] { ResultIndexOffset }, new int[] { ResultIndexOffset + Len - 1 });

                if(!(ResultIndexOffset == 0 && Len == MeanValueOT.GetLength(0)))
                    MeanValueOT = MeanValueOT.ExtractSubArrayShallow(new int[] { ResultIndexOffset }, new int[] { ResultIndexOffset + Len - 1 });
            }

            if(GradientIN != null) {
                if(GradientIN.Dimension != 3)
                    throw new ArgumentException("Dimension of gradient result array must be 3.");
                if(Len > GradientIN.GetLength(0) + ResultIndexOffset)
                    throw new ArgumentException("mismatch between Len and 0-th length of result");
                if(GradientIN.GetLength(1) != K)
                    throw new ArgumentException();
                if(GradientIN.GetLength(2) != D)
                    throw new ArgumentException();
                if(GradientOT.Dimension != 3)
                    throw new ArgumentException("Dimension of gradient result array must be 3.");
                if(Len > GradientOT.GetLength(0) + ResultIndexOffset)
                    throw new ArgumentException("mismatch between Len and 0-th length of result");
                if(GradientOT.GetLength(1) != K)
                    throw new ArgumentException();
                if(GradientOT.GetLength(2) != D)
                    throw new ArgumentException();

                if(!(ResultIndexOffset == 0 && Len == GradientIN.GetLength(0)))
                    GradientIN = GradientIN.ExtractSubArrayShallow(new int[] { ResultIndexOffset, 0, 0 }, new int[] { ResultIndexOffset + Len - 1, K - 1, D - 1 });

                if(!(ResultIndexOffset == 0 && Len == GradientOT.GetLength(0)))
                    GradientOT = GradientOT.ExtractSubArrayShallow(new int[] { ResultIndexOffset, 0, 0 }, new int[] { ResultIndexOffset + Len - 1, K - 1, D - 1 });
            }


            var coöSys = NS.GetNodeCoordinateSystem(this.GridDat);

            if(coöSys != NodeCoordinateSystem.EdgeCoord)
                throw new ArgumentOutOfRangeException();

            DGField.EvaluateEdgeInternal(e0, Len, NS, this.m_Basis, this.m_Mda_Coordinates, 
                ValueIN, ValueOT,
                MeanValueIN, MeanValueOT,
                GradientIN, GradientOT,
                ResultPreScale);
        }


        /// <summary>
        /// Evaluates the gradient of the field;
        /// </summary>
        /// <param name="j0">local index of the first cell to evaluate</param>
        /// <param name="Len">Number of cells to evaluate</param>
        /// <param name="NS">
        /// As usual, the node set.
        /// </param>
        /// <param name="result">
        /// on exit, result of the evaluations are accumulated there;
        /// the original content is scaled by <paramref name="ResultPreScale"/>;
        /// 1st index: cell index minus <paramref name="j0"/>;
        /// 2nd index: node index;
        /// 3rd index: spatial coordinate;
        /// </param>
        /// <param name="ResultCellindexOffset">
        /// an offset for the first index of <paramref name="result"/>;
        /// </param>
        /// <param name="ResultPreScale">
        /// see <paramref name="result"/>
        /// </param>
        public override void EvaluateGradient(int j0, int Len, NodeSet NS, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale) {
            int D = GridDat.SpatialDimension; // spatial dimension
            int N = m_Basis.Length;      // number of coordinates per cell
            int M = NS.NoOfNodes;        // number of nodes

            var resultAcc = result.ExtractSubArrayShallow(new int[] { ResultCellindexOffset, 0, 0 }, new int[] { ResultCellindexOffset + Len - 1, M - 1, D - 1 });
            var coord = m_Mda_Coordinates.ExtractSubArrayShallow(new int[] { j0, 0 }, new int[] { j0 + Len - 1, N - 1 });

            DGField.EvaluateGradientInternal(j0, Len, NS, this.m_Basis, coord, resultAcc, ResultPreScale);
        }

        /// <summary>
        /// Evaluates all 2nd derivatives (by cell-local analytic derivation of
        /// the basis polynomials) of this field; 
        /// </summary>
        /// <param name="j0"></param>
        /// <param name="Len"></param>
        /// <param name="NS"></param>
        /// <param name="result">
        /// <list type="bullet">
        ///   <item>1st index: cell index <em>j</em></item>
        ///   <item>2nd index: node index <em>m</em> into nodeset #<paramref name="NodeSetIndex"/></item>
        ///   <item>3rd index: spatial direction of 1st derivation, <em>k</em></item>
        ///   <item>4th index: spatial direction of 2nd derivation, <em>l</em></item>
        /// </list>
        /// So, the entry [j,m,k,l] = \f$ \frac{\partial}{\partial x_k} \frac{\partial}{\partial x_l} \varphi (\vec{\xi}_m)\f$ 
        /// where \f$ \vec{xi}_m\f$  is the <em>m</em>-th
        /// vector in the node set #<paramref name="NodeSetIndex"/>, in the
        /// <em>j</em>-th cell.
        /// </param>
        /// <remarks>
        /// Because of 2 derivatives taken, this field needs to be at least of
        /// DG degree 2 to get a non-zero result from this method.
        /// </remarks>
        public void EvaluateHessian(int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            int D = this.GridDat.SpatialDimension; // spatial dimension
            int N = m_Basis.GetLength(j0);      // number of coordinates per cell

            if (result.Dimension != 4)
                throw new ArgumentException("dimension must be 4", "result");
            if (Len > result.GetLength(0))
                throw new ArgumentOutOfRangeException("mismatch between Len and 0-th length of result");
            if (result.GetLength(2) != D || result.GetLength(3) != D)
                throw new ArgumentOutOfRangeException("3rd and 4th length of result must be " + D + "x" + D + ";");

            // not very optimized....

            MultidimensionalArray BasisHessian = this.Basis.CellEval2ndDeriv(NS,j0,Len);

            result.Multiply(1.0, this.m_Mda_Coordinates.ExtractSubArrayShallow(new int[] { j0, 0 }, new int[] { j0 + Len - 1, N - 1 }), BasisHessian, 0.0, "jkde", "jn", "jknde");


            /*
            var scales = this.GridDat.ChefBasis.Scaling;
            var Tinv = this.GridDat.Cells.InverseTransformation;

            MultidimensionalArray gbv = Basis.Evaluate2ndDeriv(NS);
            int M = gbv.GetLength(0); // number of nodes per cell
            if (result.GetLength(1) != M)
                throw new ArgumentException("length of dimension 1 must be " + M, "result");

            var Coordinates = (MultidimensionalArray)this.Coordinates;

            double[,] acc = new double[D, D];
            double[,] Tinv_j = new double[D, D];

            // loop over cells...
            for (int j = 0; j < Len; j++) {
                double sc = scales[j0 + j];

                if (this.GridDat.Cells.IsCellAffineLinear(j0 + j)) {

                    // load transformation matrix for easier access
                    for (int d1 = 0; d1 < D; d1++)
                        for (int d2 = 0; d2 < D; d2++)
                            Tinv_j[d1, d2] = Tinv[j + j0, d1, d2];

                    // loop over nodes ...
                    for (int m = 0; m < M; m++) {

                        Array.Clear(acc, 0, D * D);

                        // sum up...
                        for (int n = 0; n < N; n++) { // loop over basis functions

                            double coord = Coordinates[j + j0, n];
                            for (int d1 = 0; d1 < D; d1++) { // 1st loop over spatial dimension
                                for (int d2 = 0; d2 < D; d2++) { // 1st loop over spatial dimension
                                    acc[d1, d2] += gbv[m, n, d1, d2] * coord;
                                }
                            }
                        }

                        // transform...
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                double r = 0.0;
                                for (int k1 = 0; k1 < D; k1++) {
                                    double rr = 0;
                                    for (int k2 = 0; k2 < D; k2++) {
                                        rr += acc[k1, k2] * Tinv_j[k2, d1];
                                    }
                                    r += rr * Tinv_j[k1, d2];
                                }
                                result[j, m, d1, d2] = r * sc;
                            }
                        }
                    }
                } else {
                    throw new NotImplementedException("nonlinear cell: todo");
                }
             */
            

            /*
            
            //Test code;
             
            var HessBasis = this.Basis.CellEval2ndDeriv(NodeSetIndex, j0, Len);
            var CheckResult = result.CloneAs();
            

            for (int j = 0; j < Len; j++) {
                
                // loop over nodes ...
                for (int m = 0; m < M; m++) {

                    for (int d1 = 0; d1 < D; d1++) {
                        for (int d2 = 0; d2 < D; d2++) {
                            double accu = 0;

                            for (int n = 0; n < N; n++) {
                                accu += HessBasis[j, m, n, d1, d2]*Coordinates[j + j0, n];
                            }

                            CheckResult[j, m, d1, d2] = accu;
                        }
                    }

                }
            }

            CheckResult.Acc(-1.0, result);
            if (CheckResult.L2Norm() > 1.0e-8)
                throw new Exception();
             */
        }

        /// <summary>
        /// constructs a new field. 
        /// </summary>
        /// <param name="__Basis">The basis that is used for this field</param>
        /// <param name="__Identification">
        /// identification string for this field;
        /// This can be null or empty, 
        /// however, if IO should be preformed for this object, the identification must be unique 
        /// within the given context;
        /// </param>
        public SinglePhaseField(Basis __Basis, String __Identification)
            : base(__Basis, __Identification) {

            // allocate mem
            // ------------
            m_Mda_Coordinates = MultidimensionalArray.Create(__Basis.GridDat.iLogicalCells.Count, __Basis.MaximalLength);
        }

        /// <summary>
        /// guess what?
        /// </summary>
        /// <returns></returns>
        override public object Clone() {
            SinglePhaseField r = new SinglePhaseField(this.m_Basis);
            r.m_Identification = this.Identification; // works also if this.Identification is null!

            if (this.m_Mda_Coordinates.IsContinious && r.m_Mda_Coordinates.IsContinious) {
                Array.Copy(this.m_Mda_Coordinates.Storage, r.m_Mda_Coordinates.Storage, this.m_Mda_Coordinates.Length);
            } else {
                throw new NotImplementedException("todo");
            }

            return r;
        }

        /// <summary>
        /// see <see cref="BoSSS.Foundation.DGField.FillMPISendBuffer"/>, performance-optimized version;
        /// </summary>
        public override int FillMPISendBuffer(int proc, double[] Buffer, int st) {
            int[] CellIndexList = GridDat.iParallel.SendCommLists[proc];

            int j0 = CellIndexList[0];
            int I = CellIndexList.Length;
            int OverallLen = CellIndexList[I - 1] - CellIndexList[0] + 1;
            int N = m_Basis.MaximalLength;

            MultidimensionalArray ContStor = m_Mda_Coordinates;
            int jOffset = 0;

            if (ContStor != null && ContStor.IsContinious) {
                // accelerated version (premature opt....)
                // =======================================
                double[] storage = ContStor.Storage;

                int l = 0;
                for (int i = 0; i < I; i++) {
                    int i0 = ContStor.Index(CellIndexList[i] - jOffset, 0);
                    Array.Copy(storage, i0, Buffer, st + l, N);
                    l += N;
                }

                return l;
            } else {
                return base.FillMPISendBuffer(proc, Buffer, st);
            }
        }

        /// <summary>
        /// see <see cref="BoSSS.Foundation.DGField.CopyFromMPIrecvBuffer"/>, performance-optimized version;
        /// </summary>
        public override int CopyFromMPIrecvBuffer(int proc, double[] Buffer, int i0) {
            int j_insert = GridDat.iParallel.RcvCommListsInsertIndex[proc];
            int Len = GridDat.iParallel.RcvCommListsNoOfItems[proc];
            int N = m_Basis.MaximalLength;


            MultidimensionalArray ContStor = m_Mda_Coordinates;
            int jOffset = 0;

            if (ContStor != null && ContStor.IsContinious) {
                // accelerated version (premature opt....)
                // =======================================
                double[] storage = ContStor.Storage;
                Array.Copy(Buffer, i0, storage, ContStor.Index(j_insert - jOffset, 0), Len * N);
                return Len * N;
            } else {
                return base.CopyFromMPIrecvBuffer(proc, Buffer, i0);
            }
        }

        /*

        /// <summary>
        /// guess what?
        /// </summary>
<<<<<<< HEAD
        /// <param name="f"></param>
        /// <param name="d">
        /// 0 for the x-derivative, 1 for the y-derivative, 2 for the z-derivative
        /// </param>
        /// <param name="alpha">
        /// scaling of <paramref name="f"/>;
        /// </param>
        /// <param name="em">
        /// An optional restriction to the domain in which the derivative is computed (it may, e.g.
        /// be only required in boundary cells, so a computation over the whole domain 
        /// would be a waste of computation power. A proper execution mask for this case would be e.g. 
        /// <see cref="BoSSS.Foundation.Grid.GridData.BoundaryCells"/>.)<br/>
        /// if null, the computation is carried out in the whole domain
        /// </param>
        /// <remarks>
        /// The derivative is calculated by a cell-by-cell derivation of the DG polynomials, therefore the
        /// (effective) DG polynomial degree is one lower than the degree of <paramref name="f"/>;<br/>
        /// In comparison to <see cref="Field.DerivativeByFlux(double, Field, int, SubGrid, bool)"/>, this method should be much faster,
        /// because no quadrature is involved;
        /// </remarks>
        override public void Derivative(double alpha, Field f, int d, CellMask em) {
            using (new FuncTrace()) {
                MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);


                if (!(this.m_Basis.IsSubBasis(f.Basis) || f.Basis.IsSubBasis(this.m_Basis)))
                    throw new ArgumentException("cannot compute derivative because of incompatible basis functions.", "f");
                int D = Grid.SpatialDimension;
                if (d < 0 || d >= D)
                    throw new ArgumentException("spatial dimension out of range.", "d");

                var DGCoordinates = this.Coordinates;

                FullMatrix[,] derivMtx;
                if (this.m_Basis.MaximalLength > f.Basis.MaximalLength)
                    derivMtx = this.m_Basis.DerivativeMatrices;
                else
                    derivMtx = f.Basis.DerivativeMatrices;

                int J = GridDat.Cells.NoOfCells;
                var invTrafo = GridDat.Cells.InverseTransformation;

                FullMatrix Qd = new FullMatrix(derivMtx[0, 0].NoOfRows, derivMtx[0, 0].NoOfCols);

                IMatrix f_Coordinates = f.Coordinates;

                if (em == null) {
                    if (f.Basis.GetType() == typeof(Basis)) {
                        // other field 'f' is single phase field or some variant
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++

                        for (int j = 0; j < J; j++) {

                            int l = m_Basis.GridDat.Cells.GetRefElementIndex(j);


                            if (m_Basis.GridDat.Cells.IsCellAffineLinear(j)) {
                                // affine-linear cell
                                // ++++++++++++++++++

                                Qd.Clear();
                                for (int dd = 0; dd < D; dd++) {
                                    Qd.Acc(invTrafo[j, dd, d], derivMtx[l, dd]);
                                }
                                int N = this.Basis.GetLength(j);
                                int M = f.Basis.GetLength(j);

                                for (int n = 0; n < N; n++) {
                                    double acc = 0;
                                    for (int m = 0; m < M; m++) {
                                        acc += Qd[n, m] * f_Coordinates[j, m];
                                    }

                                    DGCoordinates[j, n] += alpha * acc;
                                }
                            } else {
                                // nonlinear/curved cell
                                // +++++++++++++++++++++
=======

         */


        /// <summary>
        /// guess what?
        /// </summary>
        new public SinglePhaseField CloneAs() {
            return (SinglePhaseField)Clone();
        }


        /// <summary>
        /// see <see cref="DGField.Acc(double,DGField,Grid.CellMask)"/>;
        /// </summary>
        public override void Acc(double mult, DGField a, Grid.CellMask cm) {
            using (new FuncTrace()) {
                if (!a.Basis.Equals(this.Basis))
                    throw new ArgumentException("Basis of 'a' must be equal to basis of this field", "a");

                SinglePhaseField _a = a as SinglePhaseField;

                if (cm == null && _a != null
                    && this.m_Mda_Coordinates.IsContinious && _a.m_Mda_Coordinates.IsContinious) {
                    // optimized branch
                    // ++++++++++++++++
                    double[] storThis = this.m_Mda_Coordinates.Storage;
                    double[] strOther = _a.m_Mda_Coordinates.Storage;

                    int i0This = this.m_Mda_Coordinates.Index(0, 0);
                    int i0Othr = _a.m_Mda_Coordinates.Index(0, 0);
                    int N = this.Mapping.LocalLength;
                    int inc = 1;

                    unsafe {
                        fixed (double* pStorThis = storThis, pStrOther = strOther) {
                            BLAS.F77_BLAS.DAXPY(ref N, ref mult, pStrOther, ref inc, pStorThis, ref inc);
                        }
                    }

                } else {
                    // default branch
                    // ++++++++++++++

                    base.Acc(mult, a, cm);
                }
            }
        }

        /// <summary>
        /// Specialized initializer for single phase fields
        /// </summary>
        [Serializable]
        public class SinglePhaseFieldInitializer : FieldInitializer {

            /// <summary>
            /// 
            /// </summary>
            /// <param name="c"></param>
            /// <returns></returns>
            public override DGField Initialize(IInitializationContext c) {
                DGField sff;
                if (c.TryGetValue(this, out sff))
                    return sff;

                var Basis = base.BasisInfo.Initialize(c);
                SinglePhaseField f = new SinglePhaseField(Basis, this.Identification);
                myInstance = f;
                c.Add(this, f);
                return f;
            }

            /// <summary>
            /// Instance of the represented object, if already created.
            /// </summary>
            [NonSerialized]
            protected SinglePhaseField myInstance;


            /// <summary>
            /// Compares the given object <paramref name="other"/> with respect
            /// to the <see cref="DGField.FieldInitializer.Identification"/>
            /// and the <see cref="DGField.FieldInitializer.BasisInfo"/>.
            /// </summary>
            public override bool Equals(Initializer<DGField> other) {
                SinglePhaseFieldInitializer initializer = other as SinglePhaseFieldInitializer;
                if (initializer == null)
                    return false;
                if (!base.BasisInfo.Equals(initializer.BasisInfo))
                    return false;
                if (!base.Identification.Equals(initializer.Identification))
                    return false;

                return true;
            }

            /// <summary>
            /// Computes a hash code based on 
            /// <see cref="DGField.FieldInitializer.Identification"/> and
            /// <see cref="DGField.FieldInitializer.BasisInfo"/>.
            /// </summary>
            public override int GetHashCode() {
                // http://stackoverflow.com/questions/1646807/quick-and-simple-hash-code-combinations
                int hash = 5; // a prime number
                hash += 73 * base.Identification.GetHashCode();
                hash += 73 * base.BasisInfo.GetHashCode();

                return hash;
            }
        }

        SinglePhaseFieldInitializer m_InstanceInfo;

        /// <summary>
        /// To support IO-architecture, NOT for direct user interaction. 
        /// Note that it is essential that this member always returns the SAME object (reference-equals)!
        /// </summary>
        public override DGField.FieldInitializer Initializer {
            get {
                if (m_InstanceInfo == null) {
                    m_InstanceInfo = new SinglePhaseFieldInitializer() {
                        BasisInfo = this.Basis.Initializer,
                        Identification = this.Identification
                    };
                }
                return m_InstanceInfo;
            }
        }
    }
}
