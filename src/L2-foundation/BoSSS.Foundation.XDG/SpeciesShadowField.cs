﻿/* =======================================================================
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
using System.Diagnostics;

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
                     owner.Identification + "#" + owner.m_CCBasis.Tracker.GetSpeciesName(specId)) {
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

            /// <summary>
            /// Setting is forbidden for shadow fields.
            /// </summary>
            public override string Identification {
                get {
                    var specId = this.m_Coordinates._SpecisId;
                    return m_Owner.Identification + "#" + m_Owner.m_CCBasis.Tracker.GetSpeciesName(specId);
                }
                set {
                    throw new NotSupportedException("Setting is forbidden -- name of shadow field is implied by owner field.");
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
            /// see <see cref="DGField.LaplacianByFlux(double,DGField,DGField,SubGrid,SubGridBoundaryModes,SubGridBoundaryModes)"/>;
            /// </summary>
            override public void LaplacianByFlux(double alpha, DGField f, DGField tmp,
                SubGrid optionalSubGrid = null,
                SubGridBoundaryModes bndMode_1stDeriv = SubGridBoundaryModes.OpenBoundary,
                SubGridBoundaryModes bndMode_2ndDeriv = SubGridBoundaryModes.OpenBoundary) {
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
                        if(ispec >= 0) { // only if species is present
                            m_Owner.Evaluate_ithSpecies(j0 + j, NodeSet, result, ResultCellindexOffset + j, ResultPreScale, ispec, ev);
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
                    throw new ArgumentOutOfRangeException("result", "dimension of result array must be 3");
                if (result.GetLength(1) != M)
                    throw new ArgumentOutOfRangeException();
                if (result.GetLength(2) != D)
                    throw new ArgumentOutOfRangeException();

                GenericEval(j0, Len, NodeSet, result, ResultCellindexOffset, ResultPreScale, M, ev);
            }

            /// <summary>
            /// Not implemented yet.
            /// </summary>
            public override void EvaluateEdge(int e0, int Len, NodeSet NS,
                MultidimensionalArray ValueIN, MultidimensionalArray ValueOT,
                MultidimensionalArray MeanValueIN, MultidimensionalArray MeanValueOT,
                MultidimensionalArray GradientIN, MultidimensionalArray GradientOT, int ResultIndexOffset, double ResultPreScale) {


                // check arguments
                // ===============
                int D = GridDat.SpatialDimension; // spatial dimension
                int N = m_Basis.Length;      // number of coordinates per cell
                int K = NS.NoOfNodes;        // number of nodes

                {
                    if ((ValueIN != null) != (ValueOT != null))
                        throw new ArgumentException();
                    if ((MeanValueIN != null) != (MeanValueOT != null))
                        throw new ArgumentException();
                    if ((GradientIN != null) != (GradientOT != null))
                        throw new ArgumentException();

                    if (ValueIN != null) {
                        if (ValueIN.Dimension != 2)
                            throw new ArgumentException("result", "dimension of result array must be 2");
                        if (Len > ValueIN.GetLength(0) + ResultIndexOffset)
                            throw new ArgumentException("mismatch between Len and 0-th length of result");
                        if (ValueIN.GetLength(1) != K)
                            throw new ArgumentException();
                        if (ValueOT.Dimension != 2)
                            throw new ArgumentException("result", "dimension of result array must be 2");
                        if (Len > ValueOT.GetLength(0) + ResultIndexOffset)
                            throw new ArgumentException("mismatch between Len and 0-th length of result");
                        if (ValueOT.GetLength(1) != K)
                            throw new ArgumentException();

                        if (!(ResultIndexOffset == 0 && Len == ValueIN.GetLength(0)))
                            ValueIN = ValueIN.ExtractSubArrayShallow(new int[] { ResultIndexOffset, 0 }, new int[] { ResultIndexOffset + Len - 1, K - 1 });

                        if (!(ResultIndexOffset == 0 && Len == ValueOT.GetLength(0)))
                            ValueOT = ValueOT.ExtractSubArrayShallow(new int[] { ResultIndexOffset, 0 }, new int[] { ResultIndexOffset + Len - 1, K - 1 });
                    }

                    if (MeanValueIN != null) {
                        if (MeanValueIN.Dimension != 1)
                            throw new ArgumentException("Dimension of mean-value result array must be 1.");
                        if (Len > MeanValueIN.GetLength(0) + ResultIndexOffset)
                            throw new ArgumentException("mismatch between Len and 0-th length of result");
                        if (MeanValueOT.Dimension != 1)
                            throw new ArgumentException("Dimension of mean-value result array must be 1.");
                        if (Len > MeanValueOT.GetLength(0) + ResultIndexOffset)
                            throw new ArgumentException("mismatch between Len and 0-th length of result");

                        if (!(ResultIndexOffset == 0 && Len == MeanValueIN.GetLength(0)))
                            MeanValueIN = MeanValueIN.ExtractSubArrayShallow(new int[] { ResultIndexOffset }, new int[] { ResultIndexOffset + Len - 1 });

                        if (!(ResultIndexOffset == 0 && Len == MeanValueOT.GetLength(0)))
                            MeanValueOT = MeanValueOT.ExtractSubArrayShallow(new int[] { ResultIndexOffset }, new int[] { ResultIndexOffset + Len - 1 });
                    }

                    if (GradientIN != null) {
                        if (GradientIN.Dimension != 3)
                            throw new ArgumentException("Dimension of gradient result array must be 3.");
                        if (Len > GradientIN.GetLength(0) + ResultIndexOffset)
                            throw new ArgumentException("mismatch between Len and 0-th length of result");
                        if (GradientIN.GetLength(1) != K)
                            throw new ArgumentException();
                        if (GradientIN.GetLength(2) != D)
                            throw new ArgumentException();
                        if (GradientOT.Dimension != 3)
                            throw new ArgumentException("Dimension of gradient result array must be 3.");
                        if (Len > GradientOT.GetLength(0) + ResultIndexOffset)
                            throw new ArgumentException("mismatch between Len and 0-th length of result");
                        if (GradientOT.GetLength(1) != K)
                            throw new ArgumentException();
                        if (GradientOT.GetLength(2) != D)
                            throw new ArgumentException();

                        if (!(ResultIndexOffset == 0 && Len == GradientIN.GetLength(0)))
                            GradientIN = GradientIN.ExtractSubArrayShallow(new int[] { ResultIndexOffset, 0, 0 }, new int[] { ResultIndexOffset + Len - 1, K - 1, D - 1 });

                        if (!(ResultIndexOffset == 0 && Len == GradientOT.GetLength(0)))
                            GradientOT = GradientOT.ExtractSubArrayShallow(new int[] { ResultIndexOffset, 0, 0 }, new int[] { ResultIndexOffset + Len - 1, K - 1, D - 1 });
                    }


                    var coöSys = NS.GetNodeCoordinateSystem(this.GridDat);

                    if (coöSys != NodeCoordinateSystem.EdgeCoord)
                        throw new ArgumentOutOfRangeException();
                }


                // real evaluation


                int[,] E2C = this.GridDat.iGeomEdges.CellIndices; // geometric edges to geometric cells
                int[][] cellL2C = this.GridDat.iLogicalCells.AggregateCellToParts; // logical cells to geometrical 
                if (cellL2C != null)
                    throw new NotImplementedException("todo");
                int[,] TrfIdx = this.GridDat.iGeomEdges.Edge2CellTrafoIndex;

                LevelSetTracker lsTrk = this.Owner.Basis.Tracker;
                var regions = lsTrk.Regions;
                var gdat = this.GridDat;
                SpeciesId mySp = this.SpeciesId;

                int i = 0; 
                while(i < Len) {

                    int ChunkLen = 1;
                    int iEdge = i + e0;

                    bool CutOrNot = IsCellCut(iEdge, E2C, regions, mySp);
                    if(CutOrNot == false) {
                        // both Neighbor cells are *not* cut
                        // standard evaluation

                        // try to find vectorization
                        while((ChunkLen + i < Len) && (IsCellCut(iEdge + ChunkLen, E2C, regions, mySp) == false)) {
                            ChunkLen++;
                        }

                        //EvaluateEdgeInternal(iEdge, ChunkLen, NS, m_Owner.Basis.NonX_Basis, 
                        //    m_Owner.m_Coordinates
                        //    )


                    } else {
                        // one or both neighbors are cut cells; 
                        // non-vectorized evaluation
                        Debug.Assert(ChunkLen == 1);
                    }                    
                    Debug.Assert(i + ChunkLen <= Len);

                    int I0 = ResultIndexOffset + i;
                    int IE = ResultIndexOffset + i + ChunkLen - 1;

                    MultidimensionalArray chunkValueIN, chunkValueOT;
                    if(ValueIN != null) {
                        chunkValueIN = ValueIN.ExtractSubArrayShallow(new[] { I0, 0 }, new[] { IE, K - 1 });
                        chunkValueOT = ValueOT.ExtractSubArrayShallow(new[] { I0, 0 }, new[] { IE, K - 1 });
                    } else {
                        chunkValueIN = ValueIN;
                        chunkValueOT = ValueOT;
                    }

                    MultidimensionalArray chunkMeanValueIN, chunkMeanValueOT;
                    if(MeanValueIN != null) {
                        chunkMeanValueIN = MeanValueIN.ExtractSubArrayShallow(new[] { I0 }, new[] { IE });
                        chunkMeanValueOT = MeanValueOT.ExtractSubArrayShallow(new[] { I0 }, new[] { IE });
                    } else {
                        chunkMeanValueIN = MeanValueIN;
                        chunkMeanValueOT = MeanValueOT;
                    }

                    MultidimensionalArray chunkGradientIN, chunkGradientOT;
                    if(GradientIN != null) {
                        chunkGradientIN = GradientIN.ExtractSubArrayShallow(new[] { I0, 0, 0 }, new[] { IE, K - 1, D - 1 });
                        chunkGradientOT = GradientOT.ExtractSubArrayShallow(new[] { I0, 0, 0 }, new[] { IE, K - 1, D - 1 });
                    } else {
                        chunkGradientIN = GradientIN;
                        chunkGradientOT = GradientOT;
                    }

                    if(CutOrNot) {
                        Debug.Assert(ChunkLen == 1);

                        int jCell0 = E2C[iEdge, 0];
                        int jCell1 = E2C[iEdge, 1];

                        NodeSet NSvol0, NSvol1;
                        NSvol0 = NS.GetVolumeNodeSet(gdat, TrfIdx[iEdge, 0], false);
                        if (jCell1 >= 0)
                            NSvol1 = NS.GetVolumeNodeSet(gdat, TrfIdx[iEdge, 1], false);
                        else
                            NSvol1 = null;
                        
                        if(ValueIN != null) {
                            this.Evaluate(jCell0, 1, NSvol0, chunkValueIN, ResultPreScale);
                        }
                        if(ValueOT != null && jCell1 >= 0) {
                            this.Evaluate(jCell1, 1, NSvol1, chunkValueOT, ResultPreScale);
                        }

                        if(MeanValueIN != null) {
                            this.EvaluateMean(jCell0, 1, chunkMeanValueIN, 0, ResultPreScale);
                        }
                        if(MeanValueOT != null && jCell1 >= 0) {
                            this.EvaluateMean(jCell1, 1, chunkMeanValueOT, 0, ResultPreScale);
                        }

                        if(GradientIN != null) {
                            this.EvaluateGradient(jCell0, 1, NSvol0, chunkGradientIN, 0, ResultPreScale);
                        }
                        if(GradientOT != null && jCell1 >= 0) {
                            this.EvaluateGradient(jCell1, 1, NSvol1, chunkGradientOT, 0, ResultPreScale);
                        }

                    } else {
                        // for edges 'iEdge' through 'iEdge + ChunkLen - 1', this shadow field is the 0-th species for IN and OT cell


                        EvaluateEdgeInternal(iEdge, ChunkLen, NS, m_Owner.Basis.NonX_Basis,
                            m_Owner.m_Coordinates.BaseStorageMda,
                            chunkValueIN, chunkValueOT,
                            chunkMeanValueIN, chunkMeanValueOT,
                            chunkGradientIN, chunkGradientOT,
                            ResultPreScale);
                    }
                    
                    i += ChunkLen;
                }
                Debug.Assert(i == Len);
            }




            static bool IsCellCut(int iEdge, int[,] E2C, LevelSetTracker.LevelSetRegions regions, SpeciesId mySp) {
                int jCell0 = E2C[iEdge, 0];

                int iSpc_cell0 = regions.GetSpeciesIndex(mySp, jCell0);
                if (iSpc_cell0 != 0)
                    return true;

                int jCell1 = E2C[iEdge, 1];
                if(jCell1 >= 0) {
                    int iSpc_cell1 = regions.GetSpeciesIndex(mySp, jCell1);
                
                    if (iSpc_cell1 != 0)
                        return true;
                }

                return false;
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

            public MatrixStructure StructureType { get; set; } = MatrixStructure.General;

            #endregion
        }
    }
}
