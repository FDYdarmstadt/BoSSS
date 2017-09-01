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
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;

namespace BoSSS.Foundation {

    /// <summary>
    /// Base class for 'normal' DG fields (<see cref="SinglePhaseField"/>,
    /// sparse fields planned, ...), i.e. fields that have a fixed number of
    /// DOF per cell.
    /// </summary>
    abstract public class ConventionalDGField : DGField {

        /// <summary>
        /// ctor
        /// </summary>
        protected ConventionalDGField(Basis b, string name) :
            base(b, name) {
            if (b.GetType() != typeof(Basis))
                throw new ArgumentException("type " + b.GetType().FullName + " is not allowed as a basis for single phase fields", "b");
        }

        /// <summary>
        /// computes the L2 - norm based on Parceval's equality
        /// </summary>
        /// <remarks>
        /// exploiting Parcevals equality for orthonormal systems, this
        /// function is (should be?) much faster than 
        /// <see cref="DGField.L2Norm(CellMask)"/>
        /// </remarks>
        public override double L2Norm(CellMask cm) {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            using (new FuncTrace()) {

                var DGCoordinates = this.Coordinates;

                double nrm2_2 = 0;
                if (cm == null) {
                    int J = GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                    int N = this.m_Basis.MaximalLength;

                    for (int j = 0; j < J; j++) {
                        for (int n = 0; n < N; n++) {
                            double k = DGCoordinates[j, n];
                            nrm2_2 += k * k;
                        }
                    }
                } else {
                    int N = this.m_Basis.MaximalLength;

                    foreach (var c in cm) {
                        int JE = c.Len + c.i0;

                        for (int j = c.i0; j < JE; j++) {
                            for (int n = 0; n < N; n++) {
                                double k = DGCoordinates[j, n];
                                nrm2_2 += k * k;
                            }
                        }
                    }
                }

                double nrmtot = 0;
                unsafe
                {
                    csMPI.Raw.Allreduce((IntPtr)(&nrm2_2), (IntPtr)(&nrmtot), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                }
                return Math.Sqrt(nrmtot);
            }
        }

        /// <summary>
        /// computes the L2 - norm based on Parceval's equality, for a specific
        /// polynomial degree <paramref name="deg"/>.
        /// </summary>
        /// <remarks>
        /// exploiting Parcevals equality for orthonormal systems, this
        /// function is (should be?) much faster than 
        /// <see cref="DGField.L2Norm(CellMask)"/>
        /// </remarks>
        public virtual double L2NormPerMode(int deg, CellMask _cm = null) {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            using (new FuncTrace()) {
                var DGCoordinates = this.Coordinates;
                IEnumerable<Chunk> cm = _cm;
                if (cm == null)
                    cm = new Chunk[] { new Chunk() { i0 = 0, Len = GridDat.iLogicalCells.NoOfLocalUpdatedCells } };

                if (this.Basis.IsOrthonormal == false)
                    throw new NotSupportedException("Not supported for non-orthonormal basis.");

                var polys = this.Basis.Polynomials;


                double nrm2_2 = 0;

                int N = this.m_Basis.MaximalLength;

                foreach (Chunk c in cm) {
                    int JE = c.Len + c.i0;

                    for (int j = c.i0; j < JE; j++) {
                        int iKref = this.GridDat.iGeomCells.GetRefElementIndex(j);

                        int N0 = polys[iKref].FirstPolynomialIndexforDegree(deg);
                        int NE = N0 + polys[iKref].NoOfPolynomialsPerDegree(deg);

                        for (int n = N0; n < NE; n++) {
                            double k = DGCoordinates[j, n];
                            nrm2_2 += k * k;
                        }
                    }
                }

                double nrmtot = 0;
                unsafe
                {
                    csMPI.Raw.Allreduce((IntPtr)(&nrm2_2), (IntPtr)(&nrmtot), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                }
                return Math.Sqrt(nrmtot);
            }
        }


        /// <summary>
        /// evaluates the mean value over a cell;
        /// of course, the mean value doesn't depend on node set or anything like that,
        /// so no information about that has to be provided.
        /// </summary>
        /// <param name="j0">local index of the first cell to evaluate</param>
        /// <param name="Len">Number of cells to evaluate</param>
        /// <param name="result">
        /// on exit, result of the evaluations are accumulated there;
        /// the original content is scaled by <paramref name="ResultPreScale"/>;
        /// 1st index: cell index minus <paramref name="j0"/>;
        /// </param>
        /// <param name="ResultCellindexOffset">
        /// an offset for the first index of <paramref name="result"/>;
        /// </param>
        /// <param name="ResultPreScale">
        /// see <paramref name="result"/>
        /// </param>
        override public void EvaluateMean(int j0, int Len, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale) {
            int D = GridDat.SpatialDimension; // spatial dimension
            if (result.Dimension != 1)
                throw new ArgumentOutOfRangeException("dimension of result array must be 1");
            if (Len > result.GetLength(0) + ResultCellindexOffset)
                throw new ArgumentOutOfRangeException("mismatch between Len and 0-th length of result");

            // we assume that the polynomial at index 0 is constant,
            // and that the mean value of all other polynomials is zero
            var Krefs = Basis.GridDat.iGeomCells.RefElements;
            MultidimensionalArray zerothOrderBasisValue = MultidimensionalArray.Create(1);
            double[] bv = new double[Krefs.Length];
            for (int l = 0; l < Krefs.Length; l++) {
                zerothOrderBasisValue.Clear();
                Debug.Assert(Basis.Polynomials[l][0].AbsoluteDegree == 0);
                Debug.Assert(Basis.Polynomials[l][0].Coeff.Length == 1);
                Basis.Polynomials[l][0].Evaluate(zerothOrderBasisValue, Krefs[l].Center);

                bv[l] = zerothOrderBasisValue[0];
            }

            MultidimensionalArray scales = Basis.Data.Scaling;
            MultidimensionalArray nonLinOrtho = null;
            for (int j = 0; j < Len; j++) {
                int iKref = Basis.GridDat.iGeomCells.GetRefElementIndex(j + j0);

                double scaling;
                if (this.GridDat.iGeomCells.IsCellAffineLinear(j + j0)) {
                    scaling = scales[j + j0];
                } else {
                    if (nonLinOrtho == null)
                        nonLinOrtho = this.GridDat.ChefBasis.OrthonormalizationTrafo.GetValue_Cell(j0, Len, 0);

                    scaling = Basis.Data.OrthonormalizationTrafo.GetValue_Cell(j, 1, 0)[0, 0, 0];
                }

                double value = 0.0;
                if (ResultPreScale != 0.0) {
                    value = result[j + ResultCellindexOffset] * ResultPreScale;
                }

                value += bv[iKref] * Coordinates[j + j0, 0] * scaling;

                result[j + ResultCellindexOffset] = value;
            }
        }

        /// <summary>
        /// see <see cref="BoSSS.Foundation.DGField.FillMPISendBuffer"/>;
        /// </summary>
        public override int FillMPISendBuffer(int proc, double[] Buffer, int st) {
            int[] CellIndexList = GridDat.iParallel.SendCommLists[proc];

            int j0 = CellIndexList[0];
            int I = CellIndexList.Length;
            int OverallLen = CellIndexList[I - 1] - CellIndexList[0] + 1;
            int N = m_Basis.MaximalLength;

            var Coordinates = this.Coordinates;

            {
                // slower version (ref Version)
                // ============================

                int l = 0;
                for (int i = 0; i < I; i++) {
                    for (int n = 0; n < N; n++)
                        Buffer[st + l + n] = Coordinates[CellIndexList[i], n];
                    l += N;
                }

                return l;
            }
        }

        /// <summary>
        /// see <see cref="BoSSS.Foundation.DGField.GetMPISendBufferSize"/>;
        /// </summary>
        public override int GetMPISendBufferSize(int p) {
            return GridDat.iParallel.SendCommLists[p].Length * m_Basis.MaximalLength;
        }

        /// <summary>
        /// see <see cref="BoSSS.Foundation.DGField.GetMPIRecvBufferSize"/>;
        /// </summary>
        public override int GetMPIRecvBufferSize(int p) {
            return GridDat.iParallel.RcvCommListsNoOfItems[p] * m_Basis.MaximalLength;
        }

        /// <summary>
        /// see <see cref="BoSSS.Foundation.DGField.CopyFromMPIrecvBuffer"/>;
        /// </summary>
        public override int CopyFromMPIrecvBuffer(int proc, double[] Buffer, int i0) {
            int j_insert = GridDat.iParallel.RcvCommListsInsertIndex[proc];
            int Len = GridDat.iParallel.RcvCommListsNoOfItems[proc];
            int N = m_Basis.MaximalLength;

            var m_Coordinates = this.Coordinates;

            {
                // slower version (ref Version)
                // ============================
                for (int j = 0; j < Len; j++) {
                    int jins = j + j_insert;

                    for (int n = 0; n < N; n++) {
                        m_Coordinates[jins, n] = Buffer[j * N + n];
                    }
                }
                return Len * N;
            }
        }


        class BrokenDerivativeForm : IVolumeForm {

            internal BrokenDerivativeForm(int d) {
                this.m_d = d;
            }

            int m_d;

            public TermActivationFlags VolTerms {
                get {
                    return TermActivationFlags.GradUxV;
                }
            }

            public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
                return GradU[0, this.m_d] * V;
            }

            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "u" };
                }
            }

            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }
        }


        /// <summary>
        /// accumulates the derivative of DG field <paramref name="f"/> 
        /// (along the <paramref name="d"/>-th axis) times <paramref name="alpha"/>
        /// to this field, i.e. <br/>
        /// this = this + <paramref name="alpha"/>* \f$ \frac{\partial}{\partial x_d} \f$ <paramref name="f"/>;
        /// </summary>
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
        override public void Derivative(double alpha, DGField f, int d, CellMask em) {

            using (var tr = new FuncTrace()) {
                MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                if (this.Basis.Degree < f.Basis.Degree - 1)
                    throw new ArgumentException("cannot compute derivative because of incompatible basis functions.", "f");
                if (f.Basis.GetType() != this.Basis.GetType())
                    throw new ArgumentException("cannot compute derivative because of incompatible basis functions.", "f");



                int D = GridDat.SpatialDimension;
                if (d < 0 || d >= D)
                    throw new ArgumentException("spatial dimension out of range.", "d");

                Quadrature.EdgeQuadratureScheme _qInsEdge;
                Quadrature.CellQuadratureScheme _qInsVol;
                {
                    _qInsEdge = (new Quadrature.EdgeQuadratureScheme(false, EdgeMask.GetEmptyMask(this.GridDat)));
                    _qInsVol = (new Quadrature.CellQuadratureScheme(true, em));
                }

                var op = (new BrokenDerivativeForm(d)).Operator(1);

                op.Evaluate(alpha, 1.0, f.Mapping, null, this.Mapping,
                    qInsEdge: _qInsEdge,
                    qInsVol: _qInsVol);

            }
        }

        /// <summary>
        /// guess what?
        /// </summary>
        new public ConventionalDGField CloneAs() {
            return (ConventionalDGField)Clone();
        }
    }
}
