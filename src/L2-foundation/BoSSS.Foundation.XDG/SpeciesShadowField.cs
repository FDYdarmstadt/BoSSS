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
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP.Utils;
using ilPSP;

namespace BoSSS.Foundation.XDG {

    public partial class XDGField {

        /// <summary>
        /// A single phase field that represents just one species (from all
        /// species in the cut-cell field); It is zero within all cells in
        /// which no DOF for the species are allocated. (The allocation is
        /// controlled by the Level set Tracker.) This field is just a shallow
        /// copy of the cut-cell-field, i.e. any manipulation on this object
        /// will be reflected onto the owning (see <see cref="Owner"/>)
        /// cut-cell-field and vice-versa.
        /// </summary>
        public class SpeciesShadowField : ConventionalDGField {

            /// <summary>
            /// ctor
            /// </summary>
            public SpeciesShadowField(XDGField owner, SpeciesId specId) :
                base(new Basis(owner.GridDat, owner.m_CCBasis.Degree),
                     owner.Identification + "-" + owner.m_CCBasis.Tracker.GetSpeciesName(specId)) {
                m_Coordinates = new SpeciesCoordinates(owner, specId);
                m_Owner = owner;
            }

            private SpeciesCoordinates m_Coordinates;

            /// <summary>
            /// DG Coordinates for the species <see cref="SpeciesId"/>.
            /// </summary>
            public override IMatrix Coordinates {
                get {
                    return m_Coordinates;
                }
            }

            private XDGField m_Owner;

            /// <summary>
            /// the linked cut-cell-field
            /// </summary>
            XDGField Owner {
                get {
                    return m_Owner;
                }
            }

            /// <summary>
            /// identification for the species that this field represents
            /// </summary>
            public SpeciesId SpeciesId {
                get {
                    return ((SpeciesCoordinates)(this.Coordinates))._SpecisId;
                }
            }

            /// <summary>
            /// the name of the species that this field represents
            /// </summary>
            public string SpeciesName {
                get {
                    return m_Owner.m_CCBasis.Tracker.GetSpeciesName(SpeciesId);
                }
            }

            /// <summary>
            ///  clone
            /// </summary>
            public override object Clone() {
                throw new NotImplementedException("it's not possible to create an non-shallow copy of a species shadow field: it refers to some other objects memory by definition.");
            }

            /// <summary>
            /// see <see cref="DGField.LaplacianByFlux(double,DGField,DGField,SubGrid,SpatialOperator.SubGridBoundaryModes,SpatialOperator.SubGridBoundaryModes)"/>;
            /// </summary>
            override public void LaplacianByFlux(double alpha, DGField f, DGField tmp,
                SubGrid optionalSubGrid = null,
                SpatialOperator.SubGridBoundaryModes bndMode_1stDeriv = SpatialOperator.SubGridBoundaryModes.OpenBoundary,
                SpatialOperator.SubGridBoundaryModes bndMode_2ndDeriv = SpatialOperator.SubGridBoundaryModes.OpenBoundary) {
                if (tmp == null)
                    // the base implementation will create the temporary field by cloning, which is a bad idea for 
                    // the species-shadow-field
                    tmp = new SinglePhaseField(this.Basis, "tmp");
                base.LaplacianByFlux(alpha, f, tmp, optionalSubGrid, bndMode_1stDeriv, bndMode_2ndDeriv);
            }

            /// <summary>
            /// see <see cref="DGField.Laplacian(double,DGField,CellMask)"/>;
            /// </summary>
            public override void Laplacian(double alpha, DGField f, CellMask em) {
                SinglePhaseField tmp = new SinglePhaseField(this.Basis, "tmp");
                int D = this.GridDat.SpatialDimension;
                for (int d = 0; d < D; d++) {
                    tmp.Clear();
                    tmp.Derivative(1.0, f, d, em);
                    this.Derivative(alpha, tmp, d, em);
                }

            }

            /// <summary>
            /// usual evaluation, just as <see cref="DGField.Evaluate(int,int,NodeSet,MultidimensionalArray,int,double)"/>;
            /// </summary>
            public override void Evaluate(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale) {

                int M = NodeSet.NoOfNodes;

                EvaluateInternalSignature ev = DGField.EvaluateInternal;


                if (result.Dimension != 2)
                    throw new ArgumentOutOfRangeException("result", "dimension of result array must be 2");
                if (result.GetLength(1) != M)
                    throw new ArgumentOutOfRangeException();

                GenericEval(j0, Len, NodeSet, result, ResultCellindexOffset, ResultPreScale, M, ev);
            }

            private void GenericEval(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result, int ResultCellindexOffset, double ResultPreScale, int M, EvaluateInternalSignature ev) {
                if (Len > result.GetLength(0) + ResultCellindexOffset)
                    throw new ArgumentOutOfRangeException("mismatch between Len and 0-th length of result");
                var lsTrk = m_Owner.Basis.Tracker;
                int j = 0;
                while (j < Len) {
                    int L = Math.Min(lsTrk.Regions.m_LenToNextChange[j0 + j], Len - j);

                    ReducedRegionCode ReducedRegionCode;
                    int NoOfSpecies;
                    NoOfSpecies = lsTrk.Regions.GetNoOfSpecies(j0 + j, out ReducedRegionCode);
                    int ispec = lsTrk.GetSpeciesIndex(ReducedRegionCode, this.SpeciesId);

                    if (NoOfSpecies == 1) {
                        if (ispec >= 0) {
                            // next "L" cells are single-phase -> vectorized evaluation
                            m_Owner.EvaluateStd(j0 + j, L, NodeSet, result, ResultCellindexOffset + j, ResultPreScale, ev);
                        } else {
                            // species is not present in given cell -> just apply prescaling to next L cells

                            int[] I0 = new int[result.Dimension];
                            I0[0] = j + ResultCellindexOffset;
                            int[] IE = result.Lengths;
                            IE[0] = I0[0] + L;
                            for (int kk = IE.Length - 1; kk >= 0; kk--)
                                IE[kk]--;

                            result.ExtractSubArrayShallow(I0, IE).Scale(ResultPreScale);
                        }
                    } else {
                        // multiphase cell
                        L = 1;
                        m_Owner.Evaluate_ithSpecies(j0 + j, NodeSet, result, ResultCellindexOffset + j, ResultPreScale, ispec, ev);
                    }
                    j += L;
                }
            }

            /// <summary>
            /// Evaluates the gradient
            /// </summary>
            public override void EvaluateGradient(int j0, int Len, NodeSet NodeSet, MultidimensionalArray result, int ResultCellindexOffset = 0, double ResultPreScale = 0.0) {
                int M = NodeSet.GetLength(0); // number of nodes per cell
                int D = GridDat.SpatialDimension;
                EvaluateInternalSignature ev = DGField.EvaluateGradientInternal;

                if (result.Dimension != 3)
                    throw new ArgumentOutOfRangeException("result", "dimension of result array must be 2");
                if (result.GetLength(1) != M)
                    throw new ArgumentOutOfRangeException();
                if (result.GetLength(2) != D)
                    throw new ArgumentOutOfRangeException();

                GenericEval(j0, Len, NodeSet, result, ResultCellindexOffset, ResultPreScale, M, ev);
            }

            /// <summary>
            /// Not implemented yet.
            /// </summary>
            public override void EvaluateEdge(int e0, int Len, NodeSet NS, MultidimensionalArray ValueIN, MultidimensionalArray ValueOT, MultidimensionalArray MeanValueIN, MultidimensionalArray MeanValueOT, MultidimensionalArray GradientIN, MultidimensionalArray GradientOT, int ResultIndexOffset, double ResultPreScale) {
                throw new NotImplementedException("todo");
            }


            /// <summary>
            /// not supported for shadow fields
            /// </summary>
            public override FieldInitializer Initializer {
                get {
                    throw new NotSupportedException("");
                }
            }

        }

        /// <summary>
        /// Used by <see cref="SpeciesShadowField"/>, a matrix that
        /// represents the DG coordinates of exactly one species;
        /// </summary>
        private class SpeciesCoordinates : IMatrix {

            /// <summary>
            /// ctor;
            /// </summary>
            internal SpeciesCoordinates(XDGField _owner, SpeciesId _SpeciesId) {
                owner = _owner;
                NSep = this.owner.m_CCBasis.DOFperSpeciesPerCell;
                LsTrk = this.owner.m_CCBasis.Tracker;
                m_SpecisId = _SpeciesId;
            }

            private XDGField owner;

            private LevelSetTracker LsTrk;

            private int NSep;

            private SpeciesId m_SpecisId;

            /// <summary>
            /// identification handle for the species id that this coordinate matrix accesses
            /// </summary>
            public SpeciesId _SpecisId {
                get {
                    return m_SpecisId;
                }
            }

            #region IMatrix Members

            /// <summary>
            /// See <see cref="IMatrix.NoOfCols"/>
            /// </summary>
            public int NoOfCols {
                get {
                    return (NSep);
                }
            }

            /// <summary>
            /// See <see cref="IMatrix.NoOfRows"/>
            /// </summary>
            public int NoOfRows {
                get {
                    return this.owner.Coordinates.NoOfRows;
                }
            }

            /// <summary>
            /// See <see cref="IMatrix"/>
            /// </summary>
            public double this[int j, int n] {
                get {
                    ReducedRegionCode rrc;
                    int NoOfSpec = LsTrk.Regions.GetNoOfSpecies(j, out rrc);

                    int iSpec = LsTrk.GetSpeciesIndex(rrc, m_SpecisId);
                    if (iSpec < 0)
                        return 0; // directly from nirvana

                    int n0;
                    n0 = iSpec * NSep;
                    if (n < NSep)
                        return owner.Coordinates[j, n + n0];
                    else
                        return owner.Coordinates[j, n + (NoOfSpec - 1) * NSep];
                }
                set {
                    ReducedRegionCode rrc;
                    int NoOfSpec = LsTrk.Regions.GetNoOfSpecies(j, out rrc);

                    int iSpec = LsTrk.GetSpeciesIndex(rrc, m_SpecisId);
                    if (iSpec < 0)
                        return; // directly to nirvana

                    int n0;
                    n0 = iSpec * NSep;
                    if (n < NSep)
                        owner.Coordinates[j, n + n0] = value;
                    else
                        owner.Coordinates[j, n + NoOfSpec * NSep] = value;
                }
            }

            /// <summary>
            /// See <see cref="IMatrix.Clear"/>
            /// </summary>
            public void Clear() {
                owner.ClearSpecies(m_SpecisId);

                // it remains to clear the external cells
                int Jext = this.owner.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                int J = this.NoOfRows;
                int N = this.NoOfCols;
                for (int j = Jext; j < J; j++) {
                    for (int n = 0; n < N; n++) {
                        this[j, n] *= 0.0;
                    }
                }
            }

            /// <summary>
            /// See <see cref="IMatrix.Scale"/>
            /// </summary>
            public void Scale(double a) {
                int Jext = this.owner.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                var ctx = owner;
                var ccBasis = owner.m_CCBasis;
                var SpecId = this.m_SpecisId;
                var coordinates = owner.Coordinates;

                int Nsep = this.NSep;

                for (int j = 0; j < Jext; j++) {
                    ReducedRegionCode rrc;
                    int NoOfSpec = ccBasis.Tracker.Regions.GetNoOfSpecies(j, out rrc); // both fields own the same tracker,
                    // so the number of species in cell j is equal for
                    // both fields, '_a' and 'this'.

                    int iSpec = ccBasis.Tracker.GetSpeciesIndex(rrc, SpecId);
                    if (iSpec < 0)
                        // species is not present in cell 'j' -> do nothing!
                        continue;

                    // species-specific DOF's:
                    int n0;
                    n0 = iSpec * Nsep;
                    for (int n = 0; n < Nsep; n++) {
                        coordinates[j, n + n0] *= a;
                    }
                }
            }

            /// <summary>
            /// accumulates <paramref name="a"/>*<paramref name="M"/> to this matrix
            /// </summary>
            public void Acc(double a, IMatrix M) {
                if (M.NoOfCols != this.NoOfCols || M.NoOfRows != this.NoOfRows)
                    throw new ArgumentException("M must have same size and length as this matrix");
                int _NoOfRows = this.NoOfRows;
                int _NoOfCols = this.NoOfCols;

                for (int i = 0; i < _NoOfRows; i++)
                    for (int j = 0; j < _NoOfCols; j++)
                        this[i, j] += a * M[i, j];
            }

            /// <summary>
            /// See <see cref="IMatrix.CopyTo(double[,],int,int)"/>
            /// </summary>
            public void CopyTo(double[,] dest, int i0, int i1) {
                int J = this.NoOfRows;
                int N = this.NoOfCols;

                for (int j = 0; j < J; j++)
                    for (int n = 0; n < N; n++)
                        dest[i0 + j, i1 + n] = this[j, n];
            }

            /// <summary>
            /// See <see cref="IMatrix.CopyTo{T}(T,bool,int)"/>
            /// </summary>
            public void CopyTo<T>(T array, bool RowWise, int arrayoffset) where T : IList<double> {
                int J = this.NoOfRows;
                int N = this.NoOfCols;

                if (RowWise) {
                    for (int j = 0; j < J; j++) {
                        for (int n = 0; n < N; n++) {
                            array[arrayoffset] = this[j, n];
                            arrayoffset++;
                        }
                    }
                } else {
                    for (int n = 0; n < N; n++) {
                        for (int j = 0; j < J; j++) {
                            array[arrayoffset] = this[j, n];
                            arrayoffset++;
                        }
                    }
                }
            }

            #endregion
        }
    }
}
