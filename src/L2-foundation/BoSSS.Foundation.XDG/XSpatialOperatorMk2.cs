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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;

using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Quadrature.FluxQuadCommon;

using static BoSSS.Foundation.SpatialOperator;


namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// An operator which is specialized in XDG Fields, i.e.
    /// it can have components which couple the phases.
    /// Mk2: enables the definition of different equation components for each phase 
    /// </summary>
    public partial class XSpatialOperatorMk2 : ISpatialOperator {

        /// <summary>
        /// A hint for implicit/nonlinear solvers, which linearization of the operator should be used
        /// </summary>
        public LinearizationHint LinearizationHint {
            get;
            set;
        }

        /// <summary>
        /// active species in this operator; all species that will try to find solutions for.
        /// </summary>
        public ICollection<string> Species {
            get {
                return m_SpeciesList.AsReadOnly();
            }
           
        }

        List<string> m_SpeciesList = new List<string>();

        private SpatialOperator FilterSpeciesOperator(ISpatialOperator op, string species, int order, EdgeQuadratureScheme eqs, CellQuadratureScheme cqs) {

            var r = new SpatialOperator(op.DomainVar, op.ParameterVar, (degDom, degParam, degCod) => order);


            foreach(string comps in op.CodomainVar) { // loop over rows
                foreach(IEquationComponent iec in op.EquationComponents[comps]) {
                    //m_SpatialOperator.EquationComponents[comps].Add(iec);
                    
                    if(iec is ISpeciesFilter fiec) {
                        string spcNmn = fiec.validSpeciesId;

                        if(!this.Species.Contains(spcNmn)) {
                            throw new ArgumentException("error in equation components for key \"" + comps + "\" SpeciesId defined in ISpeciesFilter is not given in m_Species");
                        } else {
                            r.EquationComponents[comps].Add(iec);
                        }
                    } else {
                        // no species filter defined: per default, valid for all species.
                        r.EquationComponents[comps].Add(iec);
                    }
                }
            }

            r.EdgeQuadraturSchemeProvider = g => eqs;
            r.VolumeQuadraturSchemeProvider = g => cqs;

            r.OperatorCoefficientsProvider = delegate (IGridData g, double time) {
                throw new NotImplementedException("todo");
                //this.OperatorCoefficientsProvider()
                //
                //todo();
                //return default(CoefficientSet);
            };

            r.Commit();
            return r;
        }



  

        //==========================================================================================================================
        // Taken from old XSpatialOperator
        //==========================================================================================================================
        #region XSpatialOperator

        /// <summary>
        /// Operator working on ghost edges, see e.g.
        /// @article{burman_ghost_2010,                                                             
        ///     title = {Ghost penalty},                                                        
        ///     volume = {348},                                                                 
        ///     issn = {1631073X},                                                              
        ///     url = {http://linkinghub.elsevier.com/retrieve/pii/S1631073X10002827},          
        ///     doi = {10.1016/j.crma.2010.10.006},                                             
        ///     language = {en},                                                                
        ///     number = {21-22},                                                               
        ///     urldate = {2015-09-11},                                                         
        ///     journal = {Comptes Rendus Mathematique},                                        
        ///     author = {Burman, Erik},                                                        
        ///     month = nov,                                                                    
        ///     year = {2010},                                                                  
        ///     pages = {1217--1220}
        /// }                                                
        /// </summary>
        public SpatialOperator GhostEdgesOperator {
            get;
            private set;
        }

        /// <summary>
        /// Non-coupling surface terms; originally intended to implement the flux-form of the surface tension.
        /// </summary>
        public SpatialOperator SurfaceElementOperator {
            get;
            private set;
        }

        /// <summary>
        /// edge and cell scheme for a certain species
        /// </summary>
        public struct QrSchemPair {

            /// <summary>
            /// if null, a default scheme is chosen for the integration
            /// </summary>
            public EdgeQuadratureScheme EdgeScheme;


            /// <summary>
            /// if null, a default scheme is chosen for the integration
            /// </summary>
            public CellQuadratureScheme CellScheme;
        }
        

        /// <summary>
        /// create a matrix from this operator
        /// </summary>
        public IEvaluatorLinear GetMatrixBuilder(UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) { 
            LevelSetTracker lsTrk = GetTracker(DomainVarMap.BasisS, GetBasisS(ParameterMap), CodomainVarMap.BasisS);

            return GetMatrixBuilder(lsTrk, DomainVarMap, ParameterMap, CodomainVarMap);
        }

        /// <summary>
        /// create a matrix from this operator
        /// </summary>
        public XEvaluatorLinear GetMatrixBuilder(
            LevelSetTracker lsTrk,
            UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {

            return new XEvaluatorLinear(this, lsTrk, DomainVarMap, ParameterMap, CodomainVarMap,
                1 // based on actual level-set tracker state
                );
        }

        /// <summary>
        /// explicit evaluation of the operator
        /// </summary>
        public IEvaluatorNonLin GetEvaluatorEx(IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {
            LevelSetTracker lsTrk = GetTracker(GetBasisS(DomainFields), GetBasisS(ParameterMap), CodomainVarMap.BasisS);
            return GetEvaluatorEx(lsTrk, DomainFields, ParameterMap, CodomainVarMap);
        }

        /// <summary>
        /// explicit evaluation of the operator
        /// </summary>
        public XEvaluatorNonlin GetEvaluatorEx(
            LevelSetTracker lsTrk,
            IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {

            return new XEvaluatorNonlin(this, lsTrk,
                new CoordinateMapping(DomainFields), ParameterMap, CodomainVarMap,
                1);
        }


        static Basis[] GetBasisS(IList<DGField> ParameterMap) {
            if(ParameterMap == null)
                return new Basis[0];

            return ParameterMap.Select(f => f != null ? f.Basis : default(Basis)).ToArray();
        }

        static LevelSetTracker GetTracker(IEnumerable<Basis> dom, IEnumerable<Basis> para, IEnumerable<Basis> cod) {
            LevelSetTracker lsTrk = null;
            foreach(var enu in new[] { dom, para, cod}) {
                if(enu != null) {
                    foreach(Basis b in enu) {
                        if(b != null && b is XDGBasis xb) {
                            if(lsTrk == null) {
                                lsTrk = xb.Tracker; 
                            } else {
                                if(!object.ReferenceEquals(lsTrk, xb.Tracker))
                                    throw new ArgumentException("Tracker mismatch.");
                            }
                        }
                    }
                }
            }

            if(lsTrk == null)
                throw new ArgumentException("Unable to obtain level-set tracker.");
            return lsTrk;
        }
   

        /// <summary>
        /// Computes the Jacobian matrix of the operator by finite differences.
        /// </summary>
        public IEvaluatorLinear GetFDJacobianBuilder(IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {

            LevelSetTracker lsTrk = GetTracker(GetBasisS(DomainFields), GetBasisS(ParameterMap), CodomainVarMap.BasisS);

            return GetFDJacobianBuilder(lsTrk, DomainFields, ParameterMap, CodomainVarMap);
        }

        /// <summary>
        /// Computes the Jacobian matrix of the operator by finite differences.
        /// </summary>
        public FDJacobianBuilder GetFDJacobianBuilder(
            LevelSetTracker lsTrk,
            IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) //
        {

            var xeval = this.GetEvaluatorEx(lsTrk, DomainFields, ParameterMap, CodomainVarMap);

            DelParameterUpdate delParameterUpdate = null;
            if(this.ParameterUpdate != null) {
                delParameterUpdate = this.ParameterUpdate.ParameterUpdate;
            }

            return new FDJacobianBuilder(xeval, delParameterUpdate);
        }


        /// <summary>
        /// This class acts as a frame for some other vector, and presents only those entries which are associated with a given species.
        /// </summary>
        internal class SpeciesFrameVector<V> : IList<double>
            where V : IList<double> {

            /// <summary>
            /// ctor.
            /// </summary>
            /// <param name="lsTrk_Regions">see <see cref="LevelSetTracker.Regions"/></param>
            /// <param name="spcId">species which should be framed</param>
            /// <param name="Full">the vector that should be framed</param>
            /// <param name="FullMap"></param>
            public SpeciesFrameVector(LevelSetTracker.LevelSetRegions lsTrk_Regions, SpeciesId spcId, V Full, UnsetteledCoordinateMapping FullMap)
                : this(Full, new FrameBase(lsTrk_Regions, spcId, FullMap, false)) {
            }


            /// <summary>
            /// ctor.
            /// </summary>
            public SpeciesFrameVector(V FullVec, FrameBase __fr) {
                fr = __fr;
                m_Full = FullVec;
            }

            FrameBase fr;

            V m_Full;


            /// <summary>
            /// Mapping for the framed vector.
            /// </summary>
            public UnsetteledCoordinateMapping Mapping {
                get {
                    return fr.FrameMap;
                }
            }

            /// <summary>
            ///  not supported.
            /// </summary>
            public int IndexOf(double item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// not supported.
            /// </summary>
            public void Insert(int index, double item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// not supported.
            /// </summary>
            public void RemoveAt(int index) {
                throw new NotSupportedException();
            }


            /// <summary>
            /// set/get one element
            /// </summary>
            public double this[int i] {
                get {
                    int idx = fr.Frame2Full_Loc(i);
                    return m_Full[idx];
                }
                set {
                    int idx = fr.Frame2Full_Loc(i);
                    m_Full[idx] = value;
                }
            }

            /// <summary>
            /// not supported
            /// </summary>
            public void Add(double item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// Sets all entries to 0.0.
            /// </summary>
            public void Clear() {
                int L = this.Count;
                for(int i = 0; i < L; i++)
                    this[i] = 0;
            }

            /// <summary>
            /// not supported.
            /// </summary>
            public bool Contains(double item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// copy to array.
            /// </summary>
            public void CopyTo(double[] array, int arrayIndex) {
#if DEBUG
            if(array.GetType().IsValueType)
                throw new NotSupportedException("CopyTo value type -- probably not the expected result! (Using vector struct in CopyTo(...) - operation?)");
#endif
                int L = this.Count;
                for(int i = 0; i < L; i++)
                    array[i + arrayIndex] = this[i];
            }

            /// <summary>
            /// number of elements
            /// </summary>
            public int Count {
                get {
                    return Mapping.LocalLength;
                }
            }

            /// <summary>
            /// depends on framed object
            /// </summary>
            public bool IsReadOnly {
                get {
                    return m_Full.IsReadOnly;
                }
            }

            /// <summary>
            /// not supported.
            /// </summary>
            public bool Remove(double item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// %
            /// </summary>
            public IEnumerator<double> GetEnumerator() {
                return new MyEnumerator() {
                    m_Owner = this
                };
            }

            /// <summary>
            /// %
            /// </summary>
            IEnumerator IEnumerable.GetEnumerator() {
                return new MyEnumerator() {
                    m_Owner = this
                };
            }


            class MyEnumerator : IEnumerator<double>, IEnumerator {
                internal SpeciesFrameVector<V> m_Owner;


                public double Current {
                    get {
                        int idx = m_Owner.fr.Frame2Full_Loc(cnt);
                        if(idx < 0)
                            return 0.0;
                        else
                            return m_Owner.m_Full[idx];
                    }
                }

                object IEnumerator.Current {
                    get {
                        int idx = m_Owner.fr.Frame2Full_Loc(cnt);
                        if(idx < 0)
                            return 0.0;
                        else
                            return m_Owner.m_Full[idx];
                    }
                }

                public void Dispose() {
                    m_Owner = null;
                }

                int cnt = -1;

                public bool MoveNext() {
                    cnt++;
                    return (cnt < m_Owner.Count);
                }

                public void Reset() {
                    cnt = -1;
                }
            }
        }

        /// <summary>
        /// This class acts as a frame for some other matrix, and presents only those entries which are associated with a given species.
        /// </summary>
        public class SpeciesFrameMatrix<M> : IMutableMatrixEx
            where M : IMutableMatrixEx {

            /// <summary>
            /// Pseudo-implementation.
            /// </summary>
            public object Clone() {
                throw new NotImplementedException();
            }

            /// <summary>
            /// Not Implemented.
            /// </summary>
            public void Clear() {
                throw new NotImplementedException();
            }

            /// <summary>
            /// ctor.
            /// </summary>
            /// <param name="full">the full operator matrix, from which the species <paramref name="spcId"/> should e framed (extracted)</param>
            /// <param name="lsTrk_regions">see <see cref="LevelSetTracker.Regions"/></param>
            /// <param name="spcId">the species of interest</param>
            /// <param name="fullMapRow">row mapping for the operator matrix <paramref name="full"/></param>
            /// <param name="fullMapCol">column mapping for the operator matrix <paramref name="full"/></param>
            public SpeciesFrameMatrix(M full, LevelSetTracker.LevelSetRegions lsTrk_regions, SpeciesId spcId, UnsetteledCoordinateMapping fullMapRow, UnsetteledCoordinateMapping fullMapCol)
                : this(full, new FrameBase(lsTrk_regions, spcId, fullMapRow, false), new FrameBase(lsTrk_regions, spcId, fullMapCol, true)) {
            }

            LevelSetTracker.LevelSetRegions m_LsTrk_regions;

            /// <summary>
            /// ctor.
            /// </summary>
            public SpeciesFrameMatrix(M full, FrameBase __RowFrame, FrameBase __ColFrame) {
                m_full = full;
                Debug.Assert(object.ReferenceEquals(__RowFrame.Regions, __ColFrame.Regions));
                Debug.Assert(__RowFrame.Species == __ColFrame.Species);
                RowFrame = __RowFrame;
                ColFrame = __ColFrame;
                m_LsTrk_regions = __RowFrame.Regions;

#if DEBUG
                var grdDat = RowFrame.FullMap.BasisS.First().GridDat;
                int J = grdDat.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = grdDat.iLogicalCells.Count;
                var spc = RowFrame.Species;
                Debug.Assert(RowFrame.Species.Equals(ColFrame.Species));
                var lsTrk = m_LsTrk_regions;
                Basis[] RowBase = RowMapping.BasisS.ToArray();
                Basis[] ColBase = ColMapping.BasisS.ToArray();

                var _AvailableRowIdx = new List<int>();
                for (int j = 0; j < J; j++) {
                    int iSpc = m_LsTrk_regions.GetSpeciesIndex(spc, j);

                    if (iSpc >= 0) {
                        for (int k = 0; k < RowBase.Length; k++) {
                            int N = RowBase[k].GetLength(j);
                            for (int n = 0; n < N; n++) {
                                int iRow = RowFrame.FrameMap.GlobalUniqueCoordinateIndex(k, j, n);
                                _AvailableRowIdx.Add(iRow);
                                Debug.Assert(RowFrame.Frame2Full(iRow) >= 0);
                            }
                        }
                    }
                }
                this.AvailableRowIdx = _AvailableRowIdx.ToArray();


                var _AvailableColIdx = new List<int>();
                for (int j = 0; j < JE; j++) {
                    int iSpc = m_LsTrk_regions.GetSpeciesIndex(spc, j);

                    if (iSpc >= 0) {
                        for (int k = 0; k < ColBase.Length; k++) {
                            int N = ColBase[k].GetLength(j);
                            for (int n = 0; n < N; n++) {
                                int iCol = ColFrame.FrameMap.GlobalUniqueCoordinateIndex(k, j, n);
                                _AvailableColIdx.Add(iCol);
                                Debug.Assert(ColFrame.Frame2Full(iCol) >= 0);
                            }
                        }
                    }
                }
                this.AvailableColIdx = _AvailableColIdx.ToArray();
#endif
            }

            FrameBase RowFrame;
            FrameBase ColFrame;


#if DEBUG
            int[] AvailableRowIdx;
            int[] AvailableColIdx;
#endif

            /// <summary>
            /// row mapping for the framed matrix
            /// </summary>
            public UnsetteledCoordinateMapping RowMapping {
                get {
                    return RowFrame.FrameMap;
                }
            }

            /// <summary>
            /// mpi comm
            /// </summary>
            public MPI_Comm MPI_Comm {
                get {
                    return RowPartitioning.MPI_Comm;
                }
            }

            /// <summary>
            /// column mapping for the framed matrix
            /// </summary>
            public UnsetteledCoordinateMapping ColMapping {
                get {
                    return ColFrame.FrameMap;
                }
            }

            M m_full;


            /// <summary>
            /// get a whole bunch of elements at once
            /// </summary>
            public double[] GetValues(int RowIndex, int[] ColumnIndices) {
                int L = ColumnIndices.Length;
                double[] ret = new double[L];
                for(int j = 0; j < L; j++) {
                    ret[j] = this[RowIndex, ColumnIndices[j]];
                }
                return ret;
            }

            /// <summary>
            /// set a whole bunch of elements at once
            /// </summary>
            public void SetValues(int RowIndex, int[] ColumnIndices, double[] newValues) {
                if(ColumnIndices.Length != newValues.Length)
                    throw new ArgumentException();

                int L = ColumnIndices.Length;
                for(int j = 0; j < L; j++) {
                    this[RowIndex, ColumnIndices[j]] = newValues[j];
                }
            }


            internal int iRowFrame2Full(int iFrame) {
                int iFull = RowFrame.Frame2Full(iFrame);
                return iFull;
            }

            //int iRowFull2Frame(int iFull) {
            //    return RowFrame.Full2Frame(iFull);
            //}

            internal int iColFrame2Full(int iFrame) {
                int iFull = ColFrame.Frame2Full(iFrame);
                return iFull;
            }

            //int iColFull2Frame(int iFull) {
            //    return ColFrame.Full2Frame(iFull);
            //}



            /// <summary>
            /// set/get a specific entry
            /// </summary>
            /// <param name="i">row index</param>
            /// <param name="j">column index</param>
            public double this[int i, int j] {
                get {
                    return m_full[iRowFrame2Full(i), iColFrame2Full(j)];
                }
                set {
                    m_full[iRowFrame2Full(i), iColFrame2Full(j)] = value;
                }
            }

            /// <summary>
            /// Accumulates <paramref name="Block"/>*<paramref name="alpha"/> to this matrix,
            /// at the row/column offset <paramref name="i0"/> resp. <paramref name="j0"/>.
            /// </summary>
            /// <param name="i0">Row offset.</param>
            /// <param name="j0">Column offset.</param>
            /// <param name="alpha">Scaling factor for the accumulation operation.</param>
            /// <param name="Block">Block to accumulate.</param>
            public void AccBlock(int i0, int j0, double alpha, MultidimensionalArray Block) {
                this.AccBlock(i0, j0, alpha, Block, 1.0);
            }


            /// <summary>
            /// Accumulates a block of entries to this matrix.
            /// </summary>
            /// <param name="i0">Row index offset.</param>
            /// <param name="j0">Column index offset.</param>
            /// <param name="alpha">Scaling factor for the accumulation.</param>
            /// <param name="Block">Block to add.</param>
            /// <param name="beta">pre-scaling</param>
            public void AccBlock(int i0, int j0, double alpha, MultidimensionalArray Block, double beta) {
                if(Block.Dimension != 2)
                    throw new ArgumentException();
                int I = Block.NoOfRows;
                int J = Block.NoOfCols;

                //for (int i = 0; i < I; i++)
                //    for (int j = 0; j < J; j++)
                //        this[i0 + i, j0 + j] += alpha * Block[i, j];


                if(I <= 0 || J <= 0)
                    return;

                int NoIBlk = 1;
                int[] i0S = new int[I];
                int[] i0T = new int[I];
                int[] iLT = new int[I];
                i0T[0] = iRowFrame2Full(i0);
                iLT[0] = 1;

                int NoJBlk = 1;
                int[] j0S = new int[J];
                int[] j0T = new int[J];
                int[] jLT = new int[J];
                j0T[0] = iColFrame2Full(j0);
                jLT[0] = 1;

                for(int i = 1; i < I; i++) {
                    int iT = iRowFrame2Full(i0 + i);
                    if(iT == i0T[NoIBlk - 1] + iLT[NoIBlk - 1]) {
                        iLT[NoIBlk - 1]++;
                    } else {
                        i0T[NoIBlk] = iT;
                        iLT[NoIBlk] = 1;
                        i0S[NoIBlk] = i;
                        NoIBlk++;
                    }
                }

                for(int j = 1; j < J; j++) {
                    int jT = iColFrame2Full(j0 + j);
                    if(jT == j0T[NoJBlk - 1] + jLT[NoJBlk - 1]) {
                        jLT[NoJBlk - 1]++;
                    } else {
                        j0T[NoJBlk] = jT;
                        jLT[NoJBlk] = 1;
                        j0S[NoJBlk] = j;
                        NoJBlk++;
                    }
                }


                for(int iBlk = 0; iBlk < NoIBlk; iBlk++) {
                    for(int jBlk = 0; jBlk < NoJBlk; jBlk++) {
                        var SubBlock = Block.ExtractSubArrayShallow(new int[] { i0S[iBlk], j0S[jBlk] }, new int[] { i0S[iBlk] + iLT[iBlk] - 1, j0S[jBlk] + jLT[jBlk] - 1 });
                        //double SubLinf = SubBlock.AbsSum();
                        //if(SubLinf > 0)
                        m_full.AccBlock(i0T[iBlk], j0T[jBlk], alpha, SubBlock, beta);
                    }
                }

            }

            /// <summary>
            /// depends on the framed matrix
            /// </summary>
            public bool OccupationMutable {
                get {
                    return m_full.OccupationMutable;
                }
            }

            /// <summary>
            /// read the value of the diagonal element.
            /// </summary>
            public double GetDiagonalElement(int row) {
                return this[row, row];
            }

            /// <summary>
            /// setting of diagonal element.
            /// </summary>
            public void SetDiagonalElement(int row, double val) {
                this[row, row] = val;
            }

            /// <summary>
            /// partitioning of rows over all MPI processes
            /// </summary>
            public IPartitioning RowPartitioning {
                get {
                    return RowMapping;
                }
            }

            /// <summary>
            /// partitioning of columns over all MPI processes
            /// </summary>
            public IPartitioning ColPartition {
                get {
                    return ColMapping;
                }
            }

            /// <summary>
            /// total number of rows over all MPI processes
            /// </summary>
            public int NoOfRows {
                get {
                    return (int)(RowPartitioning.TotalLength);
                }
            }

            /// <summary>
            /// total number of rows over all MPI processes
            /// </summary>
            public int NoOfCols {
                get {
                    return (int)(ColMapping.TotalLength);
                }
            }

            /// <summary>
            /// not supported
            /// </summary>
            public void SpMV<VectorType1, VectorType2>(double alpha, VectorType1 a, double beta, VectorType2 acc)
                where VectorType1 : IList<double>
                where VectorType2 : IList<double> {
                throw new NotImplementedException();
            }

            public int GetOccupiedColumnIndices(int RowIndex, ref int[] R) {
                throw new NotImplementedException();
                ////this.Ma

                //int[] FullMtxOccupied = m_full.GetOccupiedColumnIndices(iRowFrame2Full(RowIndex));

                //List<int> ret = new List<int>();
                //for (int i = 0; i < FullMtxOccupied.Length; i++) {
                //    //iColFull2Frame(FullMtxOccupied)
                //    throw new NotImplementedException();
                //}

                //return ret.ToArray();
            }

            public int GetRow(int RowIndex, ref int[] ColumnIndices, ref double[] Values) {
                //public MsrMatrix.MatrixEntry[] GetRow(int RowIndex) {
                int LR = GetOccupiedColumnIndices(RowIndex, ref ColumnIndices);
                //var row = new MsrMatrix.MatrixEntry[Occ.Length];
                if(Values == null || Values.Length < LR)
                    Values = new double[LR];

                for(int i = 0; i < LR; i++) {
                    Values[i] = this[RowIndex, ColumnIndices[i]];
                }
                return LR;
            }
        }

        #endregion


        //==========================================================================================================================
        // Reused from SpatialOperator with modifications
        //==========================================================================================================================
        #region ModSpatialOperator

        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(int,int,int,Func{int[],int[],int[],int},string[])"/>
        /// </summary>
        public XSpatialOperatorMk2(int NoOfDomFields, int NoOfParameters, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, IEnumerable<string> __Species, params string[] __varnames)
            : this(GetSubarray(__varnames, 0, NoOfDomFields), GetSubarray(__varnames, NoOfDomFields, NoOfParameters), GetSubarray(__varnames, NoOfDomFields + NoOfParameters, NoOfCodomFields), QuadOrderFunc, __Species) {
            if(NoOfCodomFields + NoOfDomFields + NoOfParameters != __varnames.Length)
                throw new ArgumentException("mismatch between number of provided variable names and given number of domain, parameter and codomain fields.");

        }

        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(int,int,Func{int[],int[],int[],int},string[])"/>
        /// </summary>
        public XSpatialOperatorMk2(int NoOfDomFields, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, IEnumerable<string> __Species, params string[] __varnames)
           : this(GetSubarray(__varnames, 0, NoOfDomFields), GetSubarray(__varnames, NoOfDomFields, NoOfCodomFields), QuadOrderFunc, __Species) {
            if(NoOfCodomFields + NoOfDomFields != __varnames.Length)
                throw new ArgumentException("mismatch between number of provided variable names and given number of domain and codomain fields.");
        }

        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(IList{string},IList{string},Func{int[],int[],int[],int})"/>
        /// </summary>
        public XSpatialOperatorMk2(IList<string> __DomainVar, IList<string> __CoDomainVar, Func<int[], int[], int[], int> QuadOrderFunc,  IEnumerable<string> __Species)
            : this(__DomainVar, null, __CoDomainVar, QuadOrderFunc, __Species) {
        }

       
        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(IList{string},IList{string},IList{string},Func{int[],int[],int[],int})"/>
        /// </summary>
        public XSpatialOperatorMk2(IList<string> __DomainVar, IList<string> __ParameterVar, IList<string> __CoDomainVar, Func<int[], int[], int[], int> QuadOrderFunc, IEnumerable<string> __Species) {
            m_DomainVar = new string[__DomainVar.Count];
            for(int i = 0; i < m_DomainVar.Length; i++) {
                if(Array.IndexOf<string>(m_DomainVar, __DomainVar[i]) >= 0)
                    throw new ArgumentException("error in domain variables list; identifier \"" + __DomainVar[i] + "\" appears twice.", "__DomainVar");
                m_DomainVar[i] = __DomainVar[i];
            }

            if(__Species == null || __Species.Count() <= 0)
                throw new ArgumentException("Empty species list.");
            if(!__Species.IsSet())
                throw new ArgumentException("Some species seems to appear more than once.");

            this.m_SpeciesList.AddRange(__Species);

            if(__ParameterVar != null) {
                m_ParameterVar = new string[__ParameterVar.Count];
                for(int i = 0; i < m_ParameterVar.Length; i++) {
                    if(Array.IndexOf<string>(m_DomainVar, __ParameterVar[i]) >= 0)
                        throw new ArgumentException("error in parameter variables list; identifier \"" + __ParameterVar[i] + "\" is already contained in the domain variables list.", "__ParameterVar");

                    if(Array.IndexOf<string>(m_ParameterVar, __ParameterVar[i]) >= 0)
                        throw new ArgumentException("error in parameter variables list; identifier \"" + __ParameterVar[i] + "\" appears twice.", "__ParameterVar");

                    m_ParameterVar[i] = __ParameterVar[i];
                }
            } else {
                m_ParameterVar = new string[0];
            }

            m_CodomainVar = new string[__CoDomainVar.Count];
            for(int i = 0; i < m_CodomainVar.Length; i++) {
                if(Array.IndexOf<string>(m_CodomainVar, __CoDomainVar[i]) >= 0)
                    throw new ArgumentException("error in codomain variables list; identifier \"" + __CoDomainVar[i] + "\" appears twice.", "__CoDomainVar");
                m_CodomainVar[i] = __CoDomainVar[i];
            }
            
            m_EquationComponents = new SortedList<string, List<IEquationComponent>>(__CoDomainVar.Count);
            foreach(var f in __CoDomainVar) {
                m_EquationComponents.Add(f, new List<IEquationComponent>());
            }
            m_EquationComponentsHelper = new _XEquationComponents(this);
            this.QuadOrderFunction = QuadOrderFunc;

            
            GhostEdgesOperator = new SpatialOperator(DomainVar, ParameterVar, CodomainVar,
                (int[] A, int[] B, int[] C) => throw new ApplicationException("should not be called - only the 'FilterSpeciesOperator(...)' should be used."));
            SurfaceElementOperator = new SpatialOperator(DomainVar, ParameterVar, CodomainVar,
                (int[] A, int[] B, int[] C) => throw new ApplicationException("should not be called - only the 'FilterSpeciesOperator(...)' should be used."));
        }

        _XEquationComponents m_EquationComponentsHelper;

        /// <summary>
        /// for each variable in <see cref="CodomainVar"/>, a
        /// collection of equation components that define the operator.
        /// 
        /// </summary>
        public IEquationComponents EquationComponents {
            get {
                return m_EquationComponentsHelper;
            }
        }


        /// <summary>
        /// implementation of <see cref="EquationComponents" />;
        /// </summary>
        public class _XEquationComponents : IEquationComponents {

            internal _XEquationComponents(XSpatialOperatorMk2 owner) {
                m_owner = owner;
            }

            XSpatialOperatorMk2 m_owner;

            /// <summary>
            /// returns the collection of equation components for one variable in the 
            /// codomain
            /// </summary>
            /// <param name="EqnName">
            /// a variable in the codomain (<see cref="SpatialOperator.CodomainVar"/>)
            /// </param>
            /// <returns></returns>
            public ICollection<IEquationComponent> this[string EqnName] {
                get {
                    if(m_owner.m_IsCommited)
                        return m_owner.m_EquationComponents[EqnName].AsReadOnly();
                    else
                        return m_owner.m_EquationComponents[EqnName];
                }
            }


            #region IEnumerable<KeyValuePair<string,IEnumerable<IEquationComponent>> Members

            /// <summary>
            /// Gets the enumerator over the equation components of the owner
            /// </summary>
            /// <returns>An enumerator</returns>
            public IEnumerator<KeyValuePair<string, IEnumerable<IEquationComponent>>> GetEnumerator() {
                return m_owner.m_EquationComponents.Select(
                    x => new KeyValuePair<string, IEnumerable<IEquationComponent>>(
                        x.Key, x.Value.AsReadOnly())).GetEnumerator();
            }

            #endregion

            #region IEnumerable Members

            /// <summary>
            /// <see cref="GetEnumerator"/>
            /// </summary>
            /// <returns>
            /// <see cref="GetEnumerator"/>
            /// </returns>
            IEnumerator IEnumerable.GetEnumerator() {
                return GetEnumerator();
            }

            #endregion
        }


        /// <summary>
        /// finalizes the assembly of the operator;
        /// Can be called only once in the lifetime of this object.
        /// After calling this method, no adding/removing of equation components is possible.
        /// </summary>
        public virtual void Commit() {
            Verify();

            if(m_IsCommited)
                throw new ApplicationException("'Commit' has already been called - it can be called only once in the lifetime of this object.");

            m_IsCommited = true;

            GhostEdgesOperator.Commit();
            SurfaceElementOperator.Commit();

        }



        #endregion

        IParameterUpdate m_ParameterUpdate;

        /// <summary>
        /// If set, used to update parameters before evaluation.
        /// </summary>
        public IParameterUpdate ParameterUpdate {
            get {
                return m_ParameterUpdate;
            }
            set {
                if(IsCommited)
                    throw new NotSupportedException("unable to change after 'Commit()'");
                m_ParameterUpdate = value;
            }
        }

        /// <summary>
        /// An operator which computes the Jacobian matrix of this operator.
        /// All components in this operator need to implement the <see cref="ISupportsJacobianComponent"/> interface in order to support this operation.
        /// </summary>
        public ISpatialOperator GetJacobiOperator(int SpatialDimension) {
            return _GetJacobiOperator(SpatialDimension);
        }


        /// <summary>
        /// An operator which computes the Jacobian matrix of this operator.
        /// All components in this operator need to implement the <see cref="ISupportsJacobianComponent"/> interface in order to support this operation.
        /// </summary>
        public XSpatialOperatorMk2 _GetJacobiOperator(int SpatialDimension) {
            if (!this.IsCommited)
                throw new InvalidOperationException("Invalid prior to calling Commit().");


            var allcomps = new List<IEquationComponent>();
            foreach (var cdo in this.CodomainVar) {
                allcomps.AddRange(this.EquationComponents[cdo]);
                allcomps.AddRange(this.GhostEdgesOperator.EquationComponents[cdo]);
                allcomps.AddRange(this.SurfaceElementOperator.EquationComponents[cdo]);
            }
            TermActivationFlags extractTaf(IEquationComponent c) {
                TermActivationFlags ret = default(TermActivationFlags);
                if (c is IVolumeForm vf) {
                    ret = ret | vf.VolTerms;
                }

                if (c is IEdgeForm ef) {
                    ret = ret | ef.BoundaryEdgeTerms;
                    ret = ret | ef.InnerEdgeTerms;
                }

                if(c is ILevelSetForm lf) {
                    ret = ret | lf.LevelSetTerms;
                }

                return ret;
            }

            var h = new JacobianParamUpdate(this.DomainVar, this.ParameterVar, allcomps, extractTaf, SpatialDimension,
                this.ParameterUpdate != null ? this.ParameterUpdate.ParameterUpdate : default(DelParameterUpdate));



            var JacobianOp = new XSpatialOperatorMk2(
                   this.DomainVar,
                   h.JacobianParameterVars,
                   this.CodomainVar,
                   this.QuadOrderFunction,
                   this.Species.ToArray());

            void CheckCoeffUpd(IEquationComponent eq, IEquationComponent eqj) {
                bool eq_suppCoeffUpd = eq is IEquationComponentCoefficient;
                bool eqj_suppCoeffUpd = eqj is IEquationComponentCoefficient;
                if (eq_suppCoeffUpd && !eqj_suppCoeffUpd)
                    throw new NotSupportedException("Form '" + eq.GetType().Name + "' supports '" + typeof(IEquationComponentCoefficient).Name + "', but Jacobian Form '" + eqj.GetType().Name + "' does not!");
            }

            foreach (string CodNmn in this.CodomainVar) {
                foreach (var eq in this.EquationComponents[CodNmn]) {
                    if (!(eq is ISupportsJacobianComponent _eq))
                        throw new NotSupportedException(string.Format("Unable to handle component {0}: To obtain a Jacobian operator, all components must implement the {1} interface.", eq.GetType().Name, typeof(ISupportsJacobianComponent).Name));
                    foreach (var eqj in _eq.GetJacobianComponents(SpatialDimension)) {
                        CheckCoeffUpd(eq, eqj);
                        JacobianOp.EquationComponents[CodNmn].Add(eqj);
                    }
                }


                foreach (var eq in this.GhostEdgesOperator.EquationComponents[CodNmn]) {
                    if (!(eq is ISupportsJacobianComponent _eq))
                        throw new NotSupportedException(string.Format("Unable to handle component {0}: To obtain a Jacobian operator, all components must implement the {1} interface.", eq.GetType().Name, typeof(ISupportsJacobianComponent).Name));
                    foreach (var eqj in _eq.GetJacobianComponents(SpatialDimension)) {
                        CheckCoeffUpd(eq, eqj);
                        JacobianOp.GhostEdgesOperator.EquationComponents[CodNmn].Add(eqj);
                    }
                }

                foreach (var eq in this.SurfaceElementOperator.EquationComponents[CodNmn]) {
                    if (!(eq is ISupportsJacobianComponent _eq))
                        throw new NotSupportedException(string.Format("Unable to handle component {0}: To obtain a Jacobian operator, all components must implement the {1} interface.", eq.GetType().Name, typeof(ISupportsJacobianComponent).Name));
                    foreach (var eqj in _eq.GetJacobianComponents(SpatialDimension)) {
                        CheckCoeffUpd(eq, eqj);
                        JacobianOp.SurfaceElementOperator.EquationComponents[CodNmn].Add(eqj);
                    }
                }

            }

            JacobianOp.ParameterUpdate = h;

            JacobianOp.Commit();
            return JacobianOp;
        }

        #region QuadSchemeProvider

        static EdgeQuadratureScheme DefaultEQSprovider(LevelSetTracker lsTrk, SpeciesId spc, XQuadSchemeHelper SchemeHelper, int quadOrder, int TrackerHistory) {
            var edgeScheme = SchemeHelper.GetEdgeQuadScheme(spc);
            return edgeScheme;
        }


        Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, EdgeQuadratureScheme> m_EdgeQuadraturSchemeProvider = DefaultEQSprovider;

        /// <summary>
        /// User-customizable factory, to specify the edge quadrature, see also <see cref="QuadOrderFunction"/>
        /// - 1st argument: current level-set tracker
        /// - 2nd argument: species which should be integrated, one of <see cref="Species"/>
        /// - 3rd argument: a default <see cref="XQuadSchemeHelper"/>
        /// - 4th argument: quadrature order
        /// - 5th argument: level-set resp. tracker history.
        /// - return: quadrature scheme        
        /// </summary>
        public Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, EdgeQuadratureScheme> EdgeQuadraturSchemeProvider {
            get {
                if(m_EdgeQuadraturSchemeProvider == null)
                    m_EdgeQuadraturSchemeProvider = DefaultEQSprovider;
                return m_EdgeQuadraturSchemeProvider;
            }
            set {
                if(IsCommited)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_EdgeQuadraturSchemeProvider = value;
            }
        }

        

        static CellQuadratureScheme DefaultCQSprovider(LevelSetTracker lsTrk, SpeciesId spc, XQuadSchemeHelper SchemeHelper, int quadOrder, int TrackerHistory) {
            var volScheme = SchemeHelper.GetVolumeQuadScheme(spc);
            return volScheme;
        }

        Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, CellQuadratureScheme> m_VolumeQuadraturSchemeProvider;

        /// <summary>
        /// User-customizable factory, to specify the cell/volume quadrature, see also <see cref="QuadOrderFunction"/>
        /// - 1st argument: current level-set tracker
        /// - 2nd argument: species which should be integrated, one of <see cref="Species"/>
        /// - 3rd argument: a default <see cref="XQuadSchemeHelper"/>
        /// - 4th argument: quadrature order
        /// - 5th argument: level-set resp. tracker history.
        /// - return: quadrature scheme
        /// </summary>
        public Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, CellQuadratureScheme> VolumeQuadraturSchemeProvider {
            get {
                if(m_VolumeQuadraturSchemeProvider == null)
                    m_VolumeQuadraturSchemeProvider = DefaultCQSprovider;
                return m_VolumeQuadraturSchemeProvider;
            }
            set {
                if(IsCommited)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_VolumeQuadraturSchemeProvider = value;
            }
        }

        #endregion

        #region Ghost_quadSchemeProvider

        static EdgeQuadratureScheme DefaultGhostEQSprovider(LevelSetTracker lsTrk, SpeciesId spc, XQuadSchemeHelper SchemeHelper, int quadOrder, int TrackerHistory) {
            var edgeScheme = SchemeHelper.GetEdgeGhostScheme(spc);
            return edgeScheme;
        }

        Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, EdgeQuadratureScheme> m_GhostEdgeQuadraturSchemeProvider;

        /// <summary>
        /// User-customizable factory, to specify the edge quadrature for the ghost-edges operator <see cref="GhostEdgesOperator"/>, see also <see cref="QuadOrderFunction"/>
        /// - 1st argument: current level-set tracker
        /// - 2nd argument: species which should be integrated, one of <see cref="Species"/>
        /// - 3rd argument: a default <see cref="XQuadSchemeHelper"/>
        /// - 4th argument: quadrature order
        /// - 5th argument: level-set resp. tracker history.
        /// - return: quadrature scheme        
        /// </summary>
        public Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, EdgeQuadratureScheme> GhostEdgeQuadraturSchemeProvider {
            get {
                if(m_GhostEdgeQuadraturSchemeProvider == null)
                    m_GhostEdgeQuadraturSchemeProvider = DefaultGhostEQSprovider;
                return m_GhostEdgeQuadraturSchemeProvider;
            }
            set {
                if(IsCommited)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_GhostEdgeQuadraturSchemeProvider = value;
            }
        }

        #endregion

        #region SurfElement_quadSchemeProvider
        static EdgeQuadratureScheme DefaultSurfElementEQSprovider(LevelSetTracker lsTrk, SpeciesId spc, XQuadSchemeHelper SchemeHelper, int quadOrder, int TrackerHistory) {
            var edgeScheme = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(spc);
            return edgeScheme;
        }


        Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, EdgeQuadratureScheme> m_SurfaceElementEdgeQuadraturSchemeProvider;

        /// <summary>
        /// User-customizable factory, to specify the edge quadrature for the <see cref="SurfaceElementOperator"/>, see also <see cref="QuadOrderFunction"/>
        /// - 1st argument: current level-set tracker
        /// - 2nd argument: species which should be integrated, one of <see cref="Species"/>
        /// - 3rd argument: a default <see cref="XQuadSchemeHelper"/>
        /// - 4th argument: quadrature order
        /// - 5th argument: level-set resp. tracker history.
        /// - return: quadrature scheme        
        /// </summary>
        public Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, EdgeQuadratureScheme> SurfaceElement_EdgeQuadraturSchemeProvider {
            get {
                if(m_SurfaceElementEdgeQuadraturSchemeProvider == null)
                    m_SurfaceElementEdgeQuadraturSchemeProvider = DefaultSurfElementEQSprovider;
                return m_SurfaceElementEdgeQuadraturSchemeProvider;
            }
            set {
                if(IsCommited)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_SurfaceElementEdgeQuadraturSchemeProvider = value;
            }
        }

        CellQuadratureScheme DefaultSurfElmCQSprovider(LevelSetTracker lsTrk, SpeciesId spc, XQuadSchemeHelper SchemeHelper, int quadOrder, int TrackerHistory) {
            var volScheme = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(spc);
            return volScheme;
        }

        Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, CellQuadratureScheme> m_SurfaceElement_VolumeQuadraturSchemeProvider;

        /// <summary>
        /// User-customizable factory, to specify the cell/volume quadrature, see also <see cref="QuadOrderFunction"/>
        /// - 1st argument: current level-set tracker
        /// - 2nd argument: species which should be integrated, one of <see cref="Species"/>
        /// - 3rd argument: a default <see cref="XQuadSchemeHelper"/>
        /// - 4th argument: quadrature order
        /// - 5th argument: level-set resp. tracker history.
        /// - return: quadrature scheme
        /// </summary>
        public Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, CellQuadratureScheme> SurfaceElement_VolumeQuadraturSchemeProvider {
            get {
                if(m_SurfaceElement_VolumeQuadraturSchemeProvider == null)
                    m_SurfaceElement_VolumeQuadraturSchemeProvider = DefaultSurfElmCQSprovider;
                return m_VolumeQuadraturSchemeProvider;
            }
            set {
                if(IsCommited)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_SurfaceElement_VolumeQuadraturSchemeProvider = value;
            }
        }
        #endregion


        DelOperatorCoefficientsProvider m_OperatorCoefficientsProvider = delegate (LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {

            var r = new CoefficientSet() {
                GrdDat = lstrk.GridDat
            };
            /*
            if(g is Grid.Classic.GridData cgdat) {
                r.CellLengthScales = cgdat.Cells.CellLengthScale;
                r.EdgeLengthScales = cgdat.Edges.h_min_Edge;

            } else {
                Console.Error.WriteLine("Rem: still missing cell length scales for grid type " + g.GetType().FullName);
            }

            todo();
            */

            throw new NotImplementedException("todo");

            //return r;
        };

        /// <summary>
        /// 
        /// </summary>
        /// <param name="lstrk">current level-set tracker</param>
        /// <param name="spc">species which should be integrated, one of <see cref="Species"/></param>
        /// <param name="quadOrder">quadrature order used for the operator (may affect length scales, etc.)</param>
        /// <param name="TrackerHistoryIdx">level-set resp. tracker history</param>
        /// <param name="time">time</param>
        /// <returns></returns>
        public delegate CoefficientSet DelOperatorCoefficientsProvider(LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time);

        /// <summary>
        /// User-customizable factory, to modify single values (e.g. Reynolds numbers)
        /// within the operator components (those implementing <see cref="IEquationComponentCoefficient"/>)
        /// Auxiliary data passed to equation components which implement <see cref="IEquationComponentCoefficient"/>.
        /// </summary>
        public DelOperatorCoefficientsProvider OperatorCoefficientsProvider {
            get {
                return m_OperatorCoefficientsProvider;
            }
            set {
                 if(IsCommited)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_OperatorCoefficientsProvider = value;
            }
        }




        //==========================================================================================================================
        // Reused from SpatialOperator without modifications
        //==========================================================================================================================
        #region SpatialOperator

        Func<int[], int[], int[], int> m_QuadOrderFunction;

        /// <summary>
        /// Function Mapping from Domain Variable Degrees, Parameter Degrees and CoDomain Variable Degrees to the Quadrature Order
        /// </summary>
        public Func<int[], int[], int[], int> QuadOrderFunction {
            get {
                return m_QuadOrderFunction;
            }
            set {
                if(IsCommited)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_QuadOrderFunction = value;
            }
        }











        static string[] GetSubarray(string[] A, int i0, int len) {
            string[] r = new string[len];
            Array.Copy(A, i0, r, 0, len);
            return r;
        }

        /// <summary>
        /// verifies all equation components;
        /// </summary>
        /// <remarks>
        /// if a component has an illegal configuration (e.g. it's arguments
        /// (<see cref="BoSSS.Foundation.IEquationComponent.ArgumentOrdering"/>) are not contained
        /// in the domain variable list (<see cref="DomainVar"/>)), an 
        /// exception is thrown;
        /// </remarks>
        internal protected void Verify() {
            foreach(var comps in m_EquationComponents.Values) {
                foreach(IEquationComponent c in comps) {
                    foreach(string varname in c.ArgumentOrdering) {
                        if(Array.IndexOf<string>(m_DomainVar, varname) < 0)
                            throw new ApplicationException("configuration error in spatial differential operator; some equation component depends on variable \""
                                + varname
                                + "\", but this name is not a member of the domain variable list.");
                    }

                    if(c.ParameterOrdering != null) {
                        foreach(string varname in c.ParameterOrdering) {
                            if(Array.IndexOf<string>(m_ParameterVar, varname) < 0)
                                throw new ApplicationException("configuration error in spatial differential operator; some equation component depends on (parameter) variable \""
                                    + varname
                                    + "\", but this name is not a member of the parameter variable list.");

                            if(c.ArgumentOrdering.Contains(varname))
                                throw new ApplicationException("configuration error in spatial differential operator; some equation component contains variable \""
                                    + varname
                                    + "\" in parameter and argument list; this is not allowed.");
                        }
                    }
                }
            }
        }


        /// <summary>
        /// Evaluation of the <see cref="QuadOrderFunction"/>.
        /// </summary>
        /// <param name="DomainMap"></param>
        /// <param name="Parameters"></param>
        /// <param name="CodomainMap"></param>
        /// <returns></returns>
        public int GetOrderFromQuadOrderFunction(IEnumerable<Basis> DomainBasis, IEnumerable<Basis> ParameterBasis, IEnumerable<Basis> CodomainBasis) {
            /// Compute Quadrature Order
            int order;
            int[] DomainDegrees = DomainBasis.Select(f => f.Degree).ToArray();
            int[] CodomainDegrees = CodomainBasis.Select(f => f.Degree).ToArray();
            int[] ParameterDegrees;
            if(ParameterBasis != null && ParameterBasis.Count() != 0) {
                ParameterDegrees = ParameterBasis.Select(b => b == null ? b.Degree : 0).ToArray();
            } else {
                ParameterDegrees = new int[] { 0 };
            };
            order = QuadOrderFunction(DomainDegrees, ParameterDegrees, CodomainDegrees);
            return order;
        }

        private static IGridData CheckArguments(UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap) {
            var GridDat = DomainMap.GridDat;
            if(!object.ReferenceEquals(GridDat, CodomainMap.GridDat))
                throw new ArgumentException("Domain and codomain map must be assigend to the same grid.");
            if(Parameters != null)
                foreach(var prm in Parameters)
                    if(prm != null && (!object.ReferenceEquals(prm.GridDat, GridDat)))
                        throw new ArgumentException(string.Format("parameter field {0} is assigned to a different grid.", prm.Identification));
            return GridDat;
        }

        /// <summary>
        /// for some codomain variable <paramref name="CodomVar"/>,
        /// this method collects 
        /// all variables in the domain (see <see cref="DomainVar"/>)
        /// on which <paramref name="CodomVar"/> depends on.
        /// </summary>
        /// <param name="CodomVar"></param>
        /// <returns>
        /// a sub-list of <see cref="DomainVar"/>;
        /// </returns>
        /// <remarks>
        /// This method invokes <see cref="Verify"/>;<br/>
        /// the returned list is in the same order as 
        /// </remarks>
        public IList<string> CollectDependentVariables(string CodomVar) {
            if(Array.IndexOf<string>(m_CodomainVar, CodomVar) < 0)
                throw new ArgumentException("the provided variable name \""
                    + CodomVar
                    + "\" is not a member of the Codomain variable list of this spatial differential operator");
            Verify();

            var comps = m_EquationComponents[CodomVar];
            List<string> ret = new List<string>();
            for(int i = 0; i < m_DomainVar.Length; i++) {
                string varName = m_DomainVar[i];

                bool isContained = false;
                foreach(var EqnComponent in comps) {
                    foreach(string name in EqnComponent.ArgumentOrdering) {
                        if(name.Equals(varName))
                            isContained = true;
                    }
                }

                if(isContained)
                    ret.Add(varName);
            }

            return ret.ToArray();
        }

        

        /// <summary>
        /// <see cref="EquationComponents"/>
        /// </summary>
        SortedList<string, List<IEquationComponent>> m_EquationComponents;


        bool m_IsCommited = false;

        /// <summary>
        /// indicates whether the equation-assembly has been finished (by calling <see cref="Commit"/>)
        /// or not.
        /// </summary>
        public bool IsCommited {
            get {
                return m_IsCommited;
            }
        }

        /// <summary>
        /// total number of equation components, in all codomain variables
        /// </summary>
        public int TotalNoOfComponents {
            get {
                return this.m_EquationComponents.Values.Sum(x => x.Count);
            }
        }



        string[] m_DomainVar;

        /// <summary>
        /// names of (DG-) variables that represent the domain of this  differential operator;
        /// These names/strings should not be confused with field identification strings
        /// (<see cref="DGField.Identification"/>), they have nothing to do with that.
        /// </summary>
        public IList<string> DomainVar {
            get {
                return (string[])m_DomainVar.Clone();
            }
        }

        string[] m_CodomainVar;

        /// <summary>
        /// names of (DG-) variables that represent the Co-Domain of this differential operator
        /// These names/strings should not be confused with field identification strings
        /// (<see cref="DGField.Identification"/>), they have nothing to do with that.
        /// </summary>
        public IList<string> CodomainVar {
            get {
                return (string[])m_CodomainVar.Clone();
            }
        }

        string[] m_ParameterVar = new string[0];

        /// <summary>
        /// names of (DG-) variables which act as parameters; 
        /// Their role is pretty similar to those of the domain variables, and for nonlinear
        /// fluxes, there is no difference.
        /// However, for <em>linear</em> fluxes, they can be used to provide some 
        /// space-depended properties as DG-fields, e.g. for providing boundary conditions
        /// or if the linear operator is some linearization of some nonlinear operator.
        /// </summary>
        public IList<string> ParameterVar {
            get {
                return (string[])m_ParameterVar.Clone();
            }
        }



        /// <summary>
        /// Returns true, if at least one of the equation components 
        /// for codomain variable <paramref name="CodomVar"/>
        /// is of some type in <paramref name="t"/>.
        /// </summary>
        public bool ContainesComponentType_PerCodomainVar(string CodomVar, params Type[] t) {
            foreach(object o in EquationComponents[CodomVar]) {
                Type[] interfaces = o.GetType().GetInterfaces();

                for(int i = 0; i < t.Length; i++)
                    if(Array.IndexOf<Type>(interfaces, t[i]) >= 0)
                        return true;
            }

            return false;
        }

        /// <summary>
        /// Returns true, if at least one of the equation components 
        /// of this operator
        /// is of some type in <paramref name="t"/>.
        /// </summary>
        public bool ContainesComponentType(params Type[] t) {
            foreach(var s in this.CodomainVar) {
                if(ContainesComponentType_PerCodomainVar(s, t))
                    return true;
            }

            return false;
        }

        #endregion
        

    }

}
