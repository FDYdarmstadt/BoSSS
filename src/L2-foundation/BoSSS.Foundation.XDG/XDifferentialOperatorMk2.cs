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

using static BoSSS.Foundation.DifferentialOperator;


namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// An operator which is specialized in XDG Fields, i.e.
    /// it can have components which couple the phases.
    /// Mk2: enables the definition of different equation components for each phase 
    /// </summary>
    public partial class XDifferentialOperatorMk2 : IDifferentialOperator {

        bool m_IsLinear;

        /// <summary>
        /// true, if the PDE defined by operator can entirely be solved by a linear solver,
        /// i.e. the Jacobian matrix, resp. the operator matrix does **not** depend on the linearization point.
        /// </summary>
        public bool IsLinear {
            get {
                return m_IsLinear;
            }
            set {
                if(IsCommitted)
                    throw new NotSupportedException("unable to change this after operator is committed.");
                m_IsLinear = value;
            }
        }

        bool m_FluxesAreNOTMultithreadSafe;


        /// <summary>
        /// Set to true, if **all** fluxes must be synchronized in multi-threaded execution.
        /// **This will come at a performance degeneration.**
        /// This is some lazy option: the default value is false,
        /// i.e., fluxes are not synchronized.
        /// <seealso cref="IMultitreadSafety"/>
        /// </summary>
        public bool FluxesAreNOTMultithreadSafe {
            get {
                return m_FluxesAreNOTMultithreadSafe;
            }
            set {
                if (IsCommitted)
                    throw new NotSupportedException("illegal to call after commit");
                m_FluxesAreNOTMultithreadSafe = value;
            }
        }

        /// <summary>
        /// <see cref="IDifferentialOperator.VectorFieldIndices"/>
        /// </summary>
        /// <remarks>
        /// Note: two ore more domain variable names of <see cref="DomainVar"/> are considered to be part of a vector field, 
        /// if these names have the same length but differ by **exactly one character**.
        /// </remarks>
        public IEnumerable<int[]> VectorFieldIndices {
            get {
                if(!this.IsCommitted)
                    throw new NotSupportedException("Operator must be committed first.");

                return PeriodicBoundaryUtils.GetVectorFieldIndices(this.DomainVar, 3);
            }
        }

        /// <summary>
        /// <see cref="IDifferentialOperator.SolverSafeguard"/>
        /// </summary>
        public SolverSafeguard SolverSafeguard {
            get;
            set;
        }

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

        public List<string> m_SpeciesList = new List<string>();

        private DifferentialOperator FilterSpeciesOperator(IDifferentialOperator op, LevelSetTracker lsTrk, string species, int order, EdgeQuadratureScheme eqs, CellQuadratureScheme cqs, int TrackerHistory, IDictionary<SpeciesId,MultidimensionalArray> CellLenScales, IDictionary<SpeciesId,MultidimensionalArray> EdgLenScales) {

            var r = new DifferentialOperator(op.DomainVar, op.ParameterVar, op.CodomainVar, (degDom, degParam, degCod) => order);
            r.FluxesAreNOTMultithreadSafe = op.FluxesAreNOTMultithreadSafe;
            r.UserDefinedValues.AddRange(this.UserDefinedValues[species]);

            foreach(string comps in op.CodomainVar) { // loop over rows
                foreach(IEquationComponent iec in op.EquationComponents[comps]) {
                    //m_SpatialOperator.EquationComponents[comps].Add(iec);
                    
                    if(iec is ISpeciesFilter fiec && fiec.ValidSpecies != null) {
                        string spcNmn = fiec.ValidSpecies;

                        if(!this.Species.Contains(spcNmn)) {
                            throw new ArgumentException("error in equation components for key \"" + comps + "\" SpeciesId defined in ISpeciesFilter is not given in m_Species");
                        } 
                        
                        if(species.Equals(fiec.ValidSpecies)){
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

                var c = this.OperatorCoefficientsProvider(lsTrk, lsTrk.GetSpeciesId(species), order, TrackerHistory, time);

                SpeciesId id = lsTrk.GetSpeciesId(species);
                CellLenScales.TryGetValue(id, out var cls);
                EdgLenScales.TryGetValue(id, out var els);


                c.CellLengthScales = cls;
                c.EdgeLengthScales = els;

                return c;
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
        public DifferentialOperator GhostEdgesOperator {
            get;
            private set;
        }

        /// <summary>
        /// Non-coupling surface terms; originally intended to implement the flux-form of the surface tension.
        /// **Note: This only considers the 0-th level-set.**
        /// </summary>
        public DifferentialOperator SurfaceElementOperator_Ls0 {
            get;
            private set;
        }

        /// <summary>
        /// Non-coupling contact-line terms.
        /// **Note: This only considers the 0-th level-set.**
        /// </summary>
        public DifferentialOperator ContactLineOperator_Ls0 {
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
            UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
            int lsTrkHistoryIndex = 1) {
            if(!IsCommitted)
                throw new NotSupportedException("Commit() (finishing operator assembly) must be called prior to evaluation.");

            return new XEvaluatorLinear(this, lsTrk, DomainVarMap, ParameterMap, CodomainVarMap, lsTrkHistoryIndex);
        }

        /// <summary>
        /// explicit evaluation of the operator
        /// </summary>
        public IEvaluatorNonLin GetEvaluatorEx(IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {
            LevelSetTracker lsTrk = GetTracker(GetBasisS(DomainFields), GetBasisS(ParameterMap), CodomainVarMap.BasisS);
            return GetEvaluatorEx(lsTrk, DomainFields, ParameterMap, CodomainVarMap);
        }

        /*
        public double[] Evaluate(IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap)
        {
            int L = CodomainVarMap.LocalLength;
            double[] ret = new double[L];
            var ev = GetEvaluatorEx(DomainFields, ParameterMap, CodomainVarMap);
            ev.Evaluate(1.0, 0.0, ret);
            return ret;
        }
        */

        /// <summary>
        /// explicit evaluation of the operator
        /// </summary>
        public XEvaluatorNonlin GetEvaluatorEx(
            LevelSetTracker lsTrk,
            IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
            int lsTrkHistoryIndex = 1) {
            if(!IsCommitted)
                throw new NotSupportedException("Commit() (finishing operator assembly) must be called prior to evaluation.");

            return new XEvaluatorNonlin(this, lsTrk,
                new CoordinateMapping(DomainFields), ParameterMap, CodomainVarMap,
                lsTrkHistoryIndex);
        }


        internal static Basis[] GetBasisS(IList<DGField> ParameterMap) {
            if(ParameterMap == null)
                return new Basis[0];

            return ParameterMap.Select(f => f != null ? f.Basis : default(Basis)).ToArray();
        }

        internal static LevelSetTracker GetTracker(IEnumerable<Basis> dom, IEnumerable<Basis> para, IEnumerable<Basis> cod) {
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
            IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
            int lsTrkHistoryIndex = 1) //

        {
            if(!IsCommitted)
                throw new NotSupportedException("Commit() (finishing operator assembly) must be called prior to evaluation.");

            var xeval = this.GetEvaluatorEx(lsTrk, DomainFields, ParameterMap, CodomainVarMap, lsTrkHistoryIndex);

            Action<double, IEnumerable<DGField>, IEnumerable<DGField>> ParamUpdate =
                delegate (double time, IEnumerable<DGField> DomF, IEnumerable<DGField> ParamF) {
                    this.InvokeParameterUpdate(time, DomF.ToArray(), ParamF.ToArray());
                };

            return new FDJacobianBuilder(xeval, ParamUpdate);
        }

        /// <summary>
        /// Clone Method
        /// </summary>
        /// <returns></returns>
        public XDifferentialOperatorMk2 CloneAs() {
            var ret = new XDifferentialOperatorMk2(this.DomainVar, this.CodomainVar, this.QuadOrderFunction, this.Species);

            foreach (string var in this.CodomainVar) {
                foreach (IEquationComponent comp in this.EquationComponents[var]) {
                    ret.EquationComponents[var].Add(comp);
                }
            }
            ret.AgglomerationThreshold= this.AgglomerationThreshold;
            ret.LinearizationHint=this.LinearizationHint;
            ret.IsLinear= this.IsLinear;
            ret.TemporalOperator= this.TemporalOperator;
            ret.FluxesAreNOTMultithreadSafe = this.FluxesAreNOTMultithreadSafe;
            return ret;
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

                var _AvailableRowIdx = new List<long>();
                for (int j = 0; j < J; j++) {
                    int iSpc = m_LsTrk_regions.GetSpeciesIndex(spc, j);

                    if (iSpc >= 0) {
                        for (int k = 0; k < RowBase.Length; k++) {
                            int N = RowBase[k].GetLength(j);
                            for (int n = 0; n < N; n++) {
                                long iRow = RowFrame.FrameMap.GlobalUniqueCoordinateIndex(k, j, n);
                                _AvailableRowIdx.Add(iRow);
                                Debug.Assert(RowFrame.Frame2Full(iRow) >= 0);
                            }
                        }
                    }
                }
                this.AvailableRowIdx = _AvailableRowIdx.ToArray();


                var _AvailableColIdx = new List<long>();
                for (int j = 0; j < JE; j++) {
                    int iSpc = m_LsTrk_regions.GetSpeciesIndex(spc, j);

                    if (iSpc >= 0) {
                        for (int k = 0; k < ColBase.Length; k++) {
                            int N = ColBase[k].GetLength(j);
                            for (int n = 0; n < N; n++) {
                                long iCol = ColFrame.FrameMap.GlobalUniqueCoordinateIndex(k, j, n);
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
            long[] AvailableRowIdx;
            long[] AvailableColIdx;
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
            public double[] GetValues(long RowIndex, long[] ColumnIndices) {
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
            public void SetValues(long RowIndex, long[] ColumnIndices, double[] newValues) {
                if(ColumnIndices.Length != newValues.Length)
                    throw new ArgumentException();

                int L = ColumnIndices.Length;
                for(int j = 0; j < L; j++) {
                    this[RowIndex, ColumnIndices[j]] = newValues[j];
                }
            }


            internal long iRowFrame2Full(long iFrame) {
                long iFull = RowFrame.Frame2Full(iFrame);
                return iFull;
            }

            //int iRowFull2Frame(int iFull) {
            //    return RowFrame.Full2Frame(iFull);
            //}

            internal long iColFrame2Full(long iFrame) {
                long iFull = ColFrame.Frame2Full(iFrame);
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
            public double this[long i, long j] {
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
            public void AccBlock(long i0, long j0, double alpha, MultidimensionalArray Block) {
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
            public void AccBlock(long i0, long j0, double alpha, MultidimensionalArray Block, double beta) {
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
                long[] i0T = new long[I];
                int[] iLT = new int[I];
                i0T[0] = iRowFrame2Full(i0);
                iLT[0] = 1;

                int NoJBlk = 1;
                int[] j0S = new int[J];
                long[] j0T = new long[J];
                int[] jLT = new int[J];
                j0T[0] = iColFrame2Full(j0);
                jLT[0] = 1;

                for(int i = 1; i < I; i++) {
                    long iT = iRowFrame2Full(i0 + i);
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
                    long jT = iColFrame2Full(j0 + j);
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
                        var SubBlock = Block.ExtractSubArrayShallow(
                            new int[] { i0S[iBlk], j0S[jBlk] }, 
                            new int[] { i0S[iBlk] + iLT[iBlk] - 1, j0S[jBlk] + jLT[jBlk] - 1 });

                        m_full.AccBlock(i0T[iBlk], j0T[jBlk], alpha, SubBlock, beta);
                    }
                }

            }

            /// <summary>
            /// Extracts a block of entries from this matrix and stores it in <paramref name="Block"/>
            /// </summary>
            /// <param name="i0">Row index offset.</param>
            /// <param name="j0">Column index offset.</param>
            /// <param name="Block"></param>
            public void ReadBlock(long i0, long j0, MultidimensionalArray Block) {
                if(Block.Dimension != 2)
                    throw new ArgumentException();
                int I = Block.NoOfRows;
                int J = Block.NoOfCols;

                for(int i = 0; i < I; i++)
                    for(int j = 0; j < J; j++)
                        Block[i, j] = this[i0 + i, j0 + j];
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
            public double GetDiagonalElement(long row) {
                return this[row, row];
            }

            /// <summary>
            /// setting of diagonal element.
            /// </summary>
            public void SetDiagonalElement(long row, double val) {
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
            public long NoOfRows {
                get {
                    return (int)(RowPartitioning.TotalLength);
                }
            }

            /// <summary>
            /// total number of rows over all MPI processes
            /// </summary>
            public long NoOfCols {
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

            public int GetOccupiedColumnIndices(long RowIndex, ref long[] R) {
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

            public int GetRow(long RowIndex, ref long[] ColumnIndices, ref double[] Values) {
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
        /// ctor, see <see cref="DifferentialOperator.DifferentialOperator(int,int,int,Func{int[],int[],int[],int},string[])"/>
        /// </summary>
        public XDifferentialOperatorMk2(int NoOfDomFields, int NoOfParameters, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, IEnumerable<string> __Species, params string[] __varnames)
            : this(GetSubarray(__varnames, 0, NoOfDomFields), GetSubarray(__varnames, NoOfDomFields, NoOfParameters), GetSubarray(__varnames, NoOfDomFields + NoOfParameters, NoOfCodomFields), QuadOrderFunc, __Species) {
            if(NoOfCodomFields + NoOfDomFields + NoOfParameters != __varnames.Length)
                throw new ArgumentException("mismatch between number of provided variable names and given number of domain, parameter and codomain fields.");

        }

        /// <summary>
        /// ctor, see <see cref="DifferentialOperator.DifferentialOperator(int,int,Func{int[],int[],int[],int},string[])"/>
        /// </summary>
        public XDifferentialOperatorMk2(int NoOfDomFields, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, IEnumerable<string> __Species, params string[] __varnames)
           : this(GetSubarray(__varnames, 0, NoOfDomFields), GetSubarray(__varnames, NoOfDomFields, NoOfCodomFields), QuadOrderFunc, __Species) {
            if(NoOfCodomFields + NoOfDomFields != __varnames.Length)
                throw new ArgumentException("mismatch between number of provided variable names and given number of domain and codomain fields.");
        }

        /// <summary>
        /// ctor
        /// </summary>
        public XDifferentialOperatorMk2(IList<string> __DomainVar, IList<string> __CoDomainVar, Func<int[], int[], int[], int> QuadOrderFunc, IEnumerable<string> __Species)
            : this(__DomainVar, null, __CoDomainVar, QuadOrderFunc, __Species) {
        }

        /// <summary>
        /// Almost empty constructor; Variable, Parameter, and Codomain/Equation names are specified by the 
        /// order in which equation components are added.
        /// </summary>
        public XDifferentialOperatorMk2(double __AgglomerationThreshold, params string[] species)
            : this(new string[0], new string[0], new string[0], QuadOrderFunc.NonLinear(2), species) {
        }


        /// <summary>
        /// ctor
        /// </summary>
        public XDifferentialOperatorMk2(IList<string> __DomainVar, IList<string> __ParameterVar, IList<string> __CoDomainVar, Func<int[], int[], int[], int> QuadOrderFunc, IEnumerable<string> __Species) {
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

            
            GhostEdgesOperator = new DifferentialOperator(DomainVar, ParameterVar, CodomainVar,
                (int[] A, int[] B, int[] C) => throw new ApplicationException("should not be called - only the 'FilterSpeciesOperator(...)' should be used."));
            SurfaceElementOperator_Ls0 = new DifferentialOperator(DomainVar, ParameterVar, CodomainVar,
                (int[] A, int[] B, int[] C) => throw new ApplicationException("should not be called - only the 'FilterSpeciesOperator(...)' should be used."));
            ContactLineOperator_Ls0 = new DifferentialOperator(DomainVar, ParameterVar, CodomainVar,
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

            internal _XEquationComponents(XDifferentialOperatorMk2 owner) {
                m_owner = owner;
            }

            XDifferentialOperatorMk2 m_owner;

            /// <summary>
            /// Returns the collection of equation components for one variable in the codomain;
            /// If the <paramref name="EqnName"/> is not known, and the operator is not committed yet (<see cref="DifferentialOperator.Commit"/>) a new 
            /// equation/codomain name is appended.
            /// </summary>
            /// <param name="EqnName">
            /// a variable in the codomain (<see cref="DifferentialOperator.CodomainVar"/>)
            /// </param>
            /// <returns></returns>
            public ICollection<IEquationComponent> this[string EqnName] {
                get {
                    if(m_owner.m_IsCommitted) {
                        return m_owner.m_EquationComponents[EqnName].AsReadOnly();
                    } else {
                        if(!m_owner.m_CodomainVar.Contains(EqnName)) {
                            m_owner.m_CodomainVar = m_owner.m_CodomainVar.Cat(EqnName);
                            m_owner.m_EquationComponents.Add(EqnName, new List<IEquationComponent>());
                        }
                        return m_owner.m_EquationComponents[EqnName];
                    }
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

        double m_AgglomerationThreshold;

        /// <summary>
        /// Cell agglomeration threshold, see <see cref="MultiphaseCellAgglomerator"/>
        /// </summary>
        public double AgglomerationThreshold {
            get {
                return m_AgglomerationThreshold;
            }
            set {
                if(IsCommitted)
                    throw new NotSupportedException("Not allowed to change the Agglomeration Threshold after operator is committed.");
                if(value < 0 || value > 1.0)
                    throw new ArgumentOutOfRangeException($"Agglomeration threshold must be between 0 and 1 (got {value}).");
                m_AgglomerationThreshold = value;
            }
        }


        /// <summary>
        /// finalizes the assembly of the operator;
        /// Can be called only once in the lifetime of this object.
        /// After calling this method, no adding/removing of equation components is possible.
        /// </summary>
        public virtual void Commit( bool allowVarAddition = true) {
             if(AgglomerationThreshold < 0 || AgglomerationThreshold > 1.0)
                    throw new ArgumentOutOfRangeException($"Agglomeration threshold must be between 0 and 1 (set as {AgglomerationThreshold}).");

            this.Verify(allowVarAddition);

            if(m_IsCommitted)
                throw new ApplicationException("'Commit' has already been called - it can be called only once in the lifetime of this object.");

            m_IsCommitted = true;

            GhostEdgesOperator.Commit();
            SurfaceElementOperator_Ls0.Commit();
            ContactLineOperator_Ls0.Commit();

            // sync the variable names of slave operators:
            // -------------------------------------------


            // this is required because we allow equations and variable names to be added _before_ Commit();
            DifferentialOperator SyncSlaveOp(DifferentialOperator slave, string slaveName) {
                if(!slave.IsCommitted)
                    throw new ApplicationException();

                foreach(var s in slave.CodomainVar)
                    if(!this.CodomainVar.Contains(s)) {
                        throw new NotSupportedException($"Found codomain variable {s} in {slaveName}, but not in main operator - not supported!");
                    }
                foreach(var s in slave.DomainVar)
                    if(!this.DomainVar.Contains(s)) {
                        throw new NotSupportedException($"Found domain variable {s} in {slaveName}, but not in main operator - not supported!");
                    }

                var R = new DifferentialOperator(this.DomainVar, this.ParameterVar, this.CodomainVar, slave.QuadOrderFunction);
                foreach(var eqname in slave.CodomainVar) {
                    foreach(var c in slave.EquationComponents[eqname])
                        R.EquationComponents[eqname].Add(c);
                }

                R.Commit();
                return R;
            }

            GhostEdgesOperator = SyncSlaveOp(GhostEdgesOperator, "GhostEdgesOperator");
            SurfaceElementOperator_Ls0 = SyncSlaveOp(SurfaceElementOperator_Ls0, "SurfaceElementOperator");
            ContactLineOperator_Ls0 = SyncSlaveOp(ContactLineOperator_Ls0, "ContactLineOperator");



            if (TemporalOperator != null) {
                TemporalOperator.Commit();
            }
        }



        #endregion


        List<DelPartialParameterUpdate> m_ParameterUpdates = new List<DelPartialParameterUpdate>();

        /// <summary>
        /// <see cref="IDifferentialOperator.ParameterUpdates"/>
        /// </summary>
        public ICollection<DelPartialParameterUpdate> ParameterUpdates {
            get {
                if (m_IsCommitted) {
                    return m_ParameterUpdates.AsReadOnly();
                } else {
                    return m_ParameterUpdates;
                }
            }
        }

        List<DelParameterFactory> m_ParameterFactories = new List<DelParameterFactory>();

        /// <summary>
        /// <see cref="IDifferentialOperator.ParameterFactories"/>
        /// </summary>
        public ICollection<DelParameterFactory> ParameterFactories {
            get {
                if (IsCommitted) {
                    return m_ParameterFactories.AsReadOnly();
                } else {
                    return m_ParameterFactories;
                }
            }
        }

        List<Action<double>> m_HomotopyUpdate = new List<Action<double>>();

        /// <summary>
        /// <see cref="IDifferentialOperator.HomotopyUpdate"/>
        /// </summary>
        public ICollection<Action<double>> HomotopyUpdate {
            get {
                if(m_IsCommitted) {
                    return m_HomotopyUpdate.AsReadOnly();
                } else {
                    return m_HomotopyUpdate;
                }
            }
        }

        double m_CurrentHomotopyValue = 1.0;

        /// <summary>
        /// Can be used by a (most likely nonlinear) solver to walk along the homotopy path.
        /// Setting (to a different value) fires all <see cref="HomotopyUpdate"/> events
        /// </summary>
        public double CurrentHomotopyValue {
            get {
                return m_CurrentHomotopyValue;
            }
            set {
                if(value < 0 || value > 1)
                    throw new ArgumentOutOfRangeException();
                double oldVal = m_CurrentHomotopyValue;
                m_CurrentHomotopyValue = value;
                if(oldVal != value) {
                    foreach(var action in m_HomotopyUpdate) {
                        action(value);
                    }
                }
            }
        }



        /// <summary>
        /// An operator which computes the Jacobian matrix of this operator.
        /// All components in this operator need to implement the <see cref="ISupportsJacobianComponent"/> interface in order to support this operation.
        /// </summary>
        public IDifferentialOperator GetJacobiOperator(int SpatialDimension) {
            return _GetJacobiOperator(SpatialDimension);
        }


        /// <summary>
        /// An operator which computes the Jacobian matrix of this operator.
        /// All components in this operator need to implement the <see cref="ISupportsJacobianComponent"/> interface in order to support this operation.
        /// </summary>
        public XDifferentialOperatorMk2 _GetJacobiOperator(int SpatialDimension) {
            if (!this.IsCommitted)
                throw new InvalidOperationException("Invalid prior to calling Commit().");


            var allcomps = new List<IEquationComponent>();
            foreach (var cdo in this.CodomainVar) {
                allcomps.AddRange(this.EquationComponents[cdo]);
                allcomps.AddRange(this.GhostEdgesOperator.EquationComponents[cdo]);
                allcomps.AddRange(this.SurfaceElementOperator_Ls0.EquationComponents[cdo]);
                allcomps.AddRange(this.ContactLineOperator_Ls0.EquationComponents[cdo]);
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

            var h = new JacobianParamUpdate(this.DomainVar, this.ParameterVar, allcomps, extractTaf, SpatialDimension);
                       
            var JacobianOp = new XDifferentialOperatorMk2(
                   this.DomainVar,
                   h.JacobianParameterVars,
                   this.CodomainVar,
                   this.QuadOrderFunction,
                   this.Species.ToArray());

            JacobianOp.FluxesAreNOTMultithreadSafe = this.FluxesAreNOTMultithreadSafe;

            if (this.TemporalOperator != null)
                JacobianOp.TemporalOperator = new TemporalOperatorContainer(JacobianOp, this.TemporalOperator);

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

                foreach (var eq in this.SurfaceElementOperator_Ls0.EquationComponents[CodNmn]) {
                    if (!(eq is ISupportsJacobianComponent _eq))
                        throw new NotSupportedException(string.Format("Unable to handle component {0}: To obtain a Jacobian operator, all components must implement the {1} interface.", eq.GetType().Name, typeof(ISupportsJacobianComponent).Name));
                    foreach (var eqj in _eq.GetJacobianComponents(SpatialDimension)) {
                        CheckCoeffUpd(eq, eqj);
                        JacobianOp.SurfaceElementOperator_Ls0.EquationComponents[CodNmn].Add(eqj);
                    }
                }

                foreach (var eq in this.ContactLineOperator_Ls0.EquationComponents[CodNmn]) {
                    if (!(eq is ISupportsJacobianComponent _eq))
                        throw new NotSupportedException(string.Format("Unable to handle component {0}: To obtain a Jacobian operator, all components must implement the {1} interface.", eq.GetType().Name, typeof(ISupportsJacobianComponent).Name));
                    foreach (var eqj in _eq.GetJacobianComponents(SpatialDimension)) {
                        CheckCoeffUpd(eq, eqj);
                        JacobianOp.ContactLineOperator_Ls0.EquationComponents[CodNmn].Add(eqj);
                    }
                }

            }

            foreach(string domName in this.DomainVar)
                JacobianOp.FreeMeanValue[domName] = this.FreeMeanValue[domName];
            JacobianOp.EdgeQuadraturSchemeProvider = this.EdgeQuadraturSchemeProvider;
            JacobianOp.VolumeQuadraturSchemeProvider = this.VolumeQuadraturSchemeProvider;
            JacobianOp.SurfaceElement_VolumeQuadraturSchemeProvider = this.SurfaceElement_VolumeQuadraturSchemeProvider;
            JacobianOp.SurfaceElement_EdgeQuadraturSchemeProvider = this.SurfaceElement_EdgeQuadraturSchemeProvider;
            JacobianOp.ContactLine_VolumeQuadratureSchemeProvider = this.ContactLine_VolumeQuadratureSchemeProvider;
            JacobianOp.GhostEdgeQuadraturSchemeProvider = this.GhostEdgeQuadraturSchemeProvider;
            foreach(var species in this.Species) {
                var src = this.UserDefinedValues[species];
                var dst = JacobianOp.UserDefinedValues[species];
                foreach(var kv in src)
                    dst.Add(kv);
            }
            JacobianOp.LinearizationHint = LinearizationHint.AdHoc;


            foreach (DelParameterFactory f in this.ParameterFactories)
                JacobianOp.ParameterFactories.Add(f);
            foreach (DelPartialParameterUpdate f in this.ParameterUpdates) {
                JacobianOp.ParameterUpdates.Add(f);
            }
            JacobianOp.ParameterFactories.Add(h.AllocateParameters);
            JacobianOp.ParameterUpdates.Add(h.PerformUpdate);

            JacobianOp.OperatorCoefficientsProvider = this.OperatorCoefficientsProvider;
            JacobianOp.m_HomotopyUpdate.AddRange(this.m_HomotopyUpdate);
            JacobianOp.m_CurrentHomotopyValue = this.m_CurrentHomotopyValue;

            JacobianOp.Commit(false);
            return JacobianOp;
        }

        /// <summary>
        /// Used by <see cref="_GetJacobiOperator(int)"/> to encalsulate the temporal operator
        /// of this operator (because of the ownership, the temporal operator cannot be reused).
        /// </summary>
        class TemporalOperatorContainer : ITemporalOperator {

            XDifferentialOperatorMk2 m_newOwner;
            ITemporalOperator m_encapsulatedObj;
            public TemporalOperatorContainer(XDifferentialOperatorMk2 __newOwner, ITemporalOperator __encapsulatedObj) {
                m_encapsulatedObj = __encapsulatedObj;
                m_newOwner = __newOwner;


            }

            bool m_IsCommited;

            /// <summary>
            /// locks the configuration of the operator
            /// </summary>
            public void Commit() {
                if (m_IsCommited)
                    throw new ApplicationException("'Commit' has already been called - it can be called only once in the lifetime of this object.");
                m_IsCommited = true;

            }

            public IEvaluatorLinear GetMassMatrixBuilder(UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) {
                return m_encapsulatedObj.GetMassMatrixBuilder(DomainVarMap, ParameterMap, CodomainVarMap);
            }
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
                if(IsCommitted)
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
                if(IsCommitted)
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
                if(IsCommitted)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_GhostEdgeQuadraturSchemeProvider = value;
            }
        }

        #endregion

        #region SurfElement_quadSchemeProvider
        static EdgeQuadratureScheme DefaultSurfElementEQSprovider(LevelSetTracker lsTrk, SpeciesId spc, XQuadSchemeHelper SchemeHelper, int quadOrder, int TrackerHistory) {
            var edgeScheme = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(spc, 0);
            return edgeScheme;
        }


        Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, EdgeQuadratureScheme> m_SurfaceElementEdgeQuadraturSchemeProvider;

        /// <summary>
        /// User-customizable factory, to specify the edge quadrature for the <see cref="SurfaceElementOperator_Ls0"/>, see also <see cref="QuadOrderFunction"/>
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
                if(IsCommitted)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_SurfaceElementEdgeQuadraturSchemeProvider = value;
            }
        }

        CellQuadratureScheme DefaultSurfElmCQSprovider(LevelSetTracker lsTrk, SpeciesId spc, XQuadSchemeHelper SchemeHelper, int quadOrder, int TrackerHistory) {
            var volScheme = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(spc, 0);
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
                return m_SurfaceElement_VolumeQuadraturSchemeProvider;
            }
            set {
                if(IsCommitted)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_SurfaceElement_VolumeQuadraturSchemeProvider = value;
            }
        }
        #endregion

        Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, CellQuadratureScheme> m_ContactLine_VolumeQuadraturSchemeProvider;

        /// <summary>
        /// User-customizable factory, to specify the cell/volume quadrature, see also <see cref="QuadOrderFunction"/>
        /// - 1st argument: current level-set tracker
        /// - 2nd argument: species which should be integrated, one of <see cref="Species"/>
        /// - 3rd argument: a default <see cref="XQuadSchemeHelper"/>
        /// - 4th argument: quadrature order
        /// - 5th argument: level-set resp. tracker history.
        /// - return: quadrature scheme
        /// </summary>
        public Func<LevelSetTracker, SpeciesId, XQuadSchemeHelper, int, int, CellQuadratureScheme> ContactLine_VolumeQuadratureSchemeProvider {
            get {
                if (m_ContactLine_VolumeQuadraturSchemeProvider == null)
                    m_ContactLine_VolumeQuadraturSchemeProvider = DefaultContactLineCQSprovider;
                return m_ContactLine_VolumeQuadraturSchemeProvider;
            }
            set {
                if (IsCommitted)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_ContactLine_VolumeQuadraturSchemeProvider = value;
            }
        }
        
        CellQuadratureScheme DefaultContactLineCQSprovider(LevelSetTracker lsTrk, SpeciesId spc, XQuadSchemeHelper SchemeHelper, int quadOrder, int TrackerHistory) {
            var volScheme = SchemeHelper.GetContactLineQuadScheme(spc, 0);
            return volScheme;
        }


        DelOperatorCoefficientsProvider m_OperatorCoefficientsProvider; 
        
        
        CoefficientSet DefaultOperatorCoefficientsProvider (LevelSetTracker lstrk, SpeciesId spc, int quadOrder, int TrackerHistoryIdx, double time) {

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

            foreach(var kv in UserDefinedValues[lstrk.GetSpeciesName(spc)]) {
                r.UserDefinedValues.Add(kv.Key, kv.Value);
            }

            r.HomotopyValue = this.m_CurrentHomotopyValue;

            return r;
        }

        /// <summary>
        /// internal access for hack in <see cref="DependentXTemporalOperator"/>.
        /// </summary>
        internal Dictionary<string, IDictionary<string, object>> m_UserDefinedValues;

        /// <summary>
        /// Modification of <see cref="CoefficientSet.UserDefinedValues"/>, **but only if** default setting for <see cref="OperatorCoefficientsProvider"/> is used
        /// - 1st key: species
        /// - 2nd key: user-defined name 
        /// - value: object to pass to <see cref="IEquationComponentCoefficient.CoefficientUpdate"/>
        /// </summary>
        public IReadOnlyDictionary<string, IDictionary<string, object>> UserDefinedValues {
            get {
                if(m_UserDefinedValues == null) {
                    m_UserDefinedValues = new Dictionary<string, IDictionary<string, object>>();
                    foreach(string spcName in this.Species) {
                        m_UserDefinedValues.Add(spcName, new Dictionary<string, object>());
                    }
                }
                return m_UserDefinedValues;
            }
            
        }



        /// <summary>
        /// <see cref="OperatorCoefficientsProvider"/>
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
                if(m_OperatorCoefficientsProvider == null)
                    m_OperatorCoefficientsProvider = DefaultOperatorCoefficientsProvider;
                return m_OperatorCoefficientsProvider;
            }
            set {
                //if(IsCommited)
                //    throw new NotSupportedException("not allowed to change after Commit");
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
                if(IsCommitted)
                    throw new NotSupportedException("not allowed to change after Commit");
                m_QuadOrderFunction = value;
            }
        }



      

        /// <summary>
        /// <see cref="IDifferentialOperator.IsValidDomainDegreeCombination"/>
        /// </summary>
        public bool IsValidDomainDegreeCombination(int[] DomainDegreesPerVariable, int[] CodomainDegreesPerVariable) {
            if (!this.IsCommitted)
                throw new InvalidOperationException("Invalid prior to calling Commit().");

            if (DomainDegreesPerVariable.Length != this.DomainVar.Count)
                throw new ArgumentException("Mismatch between length of input and number of domain variables.");
            if (CodomainDegreesPerVariable.Length != this.CodomainVar.Count)
                throw new ArgumentException("Mismatch between length of input and number of codomain variables.");

            int i = 0;
            foreach (var cod in this.CodomainVar) {
                foreach (var comp in this.EquationComponents[cod]) {
                    if (comp is IDGdegreeConstraint dgconstr) {
                        int[] argDegrees;
                        if (comp.ArgumentOrdering != null)
                            argDegrees = comp.ArgumentOrdering.Select(argName => DomainDegreesPerVariable[this.DomainVar.IndexWhere(domName => domName == argName)]).ToArray();
                        else
                            argDegrees = new int[0];

                        if (dgconstr.IsValidDomainDegreeCombination(argDegrees, CodomainDegreesPerVariable[i]) == false)
                            return false;
                    }
                }

                i++;
            }


            return true;
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
        /// exception is thrown, if <paramref name="allowVarAddition"/> is set true, see also <see cref="IDifferentialOperator.Commit(bool)"/>
        /// </remarks>
        internal protected void Verify(bool allowVarAddition) {
            //if(this.IsLinear && LinearizationHint != LinearizationHint.AdHoc)
            //    throw new NotSupportedException("Configuration Error: for a supposedly linear operator, the linearization hint must be " + LinearizationHint.AdHoc);

            foreach (var comps in m_EquationComponents.Values) {
                foreach (IEquationComponent c in comps) {
                    foreach (string varname in c.ArgumentOrdering) {
                        if (Array.IndexOf<string>(m_DomainVar, varname) < 0) {

                            if (allowVarAddition)
                                m_DomainVar = m_DomainVar.Cat(varname);
                            else {
                                throw new ApplicationException("configuration error in spatial differential operator; equation component " + c.ToString() + " depends on variable \""
                                    + varname
                                    + "\", but this name is not a member of the domain variable list: "
                                    + m_DomainVar.ToConcatString("[", ", ", "]"));
                            }

                        }
                    }

                    if (c.ParameterOrdering != null) {
                        foreach (string varname in c.ParameterOrdering) {
                            if (Array.IndexOf<string>(m_ParameterVar, varname) < 0) {
                                if (allowVarAddition)
                                    m_ParameterVar = m_ParameterVar.Cat(varname);
                                else {
                                    throw new ApplicationException("configuration error in spatial differential operator; equation component " + c.ToString() + " depends on (parameter) variable \""
                                        + varname
                                        + "\", but this name is not a member of the parameter variable list: "
                                        + m_ParameterVar.ToConcatString("[", ", ", "]"));
                                }
                            }
                        }
                    }
                }
            }


            foreach (var comps in m_EquationComponents.Values) {
                foreach (IEquationComponent c in comps) {
                    if (c.ParameterOrdering != null) {
                        foreach (string varname in c.ParameterOrdering) {


                            if (this.m_DomainVar.Contains(varname))
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
        public int GetOrderFromQuadOrderFunction(IEnumerable<Basis> DomainBasis, IEnumerable<Basis> ParameterBasis, IEnumerable<Basis> CodomainBasis) {
            // Compute Quadrature Order
            int order;
            int[] DomainDegrees = DomainBasis.Select(f => f.Degree).ToArray();
            int[] CodomainDegrees = CodomainBasis.Select(f => f.Degree).ToArray();
            int[] ParameterDegrees;
            if(ParameterBasis != null && ParameterBasis.Count() != 0) {
                ParameterDegrees = ParameterBasis.Select(b => b != null ? b.Degree : 0).ToArray();
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
            Verify(false);

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


        bool m_IsCommitted = false;

        /// <summary>
        /// indicates whether the equation-assembly has been finished (by calling <see cref="Commit"/>)
        /// or not.
        /// </summary>
        public bool IsCommitted {
            get {
                return m_IsCommitted;
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
                return m_DomainVar.ToList().AsReadOnly();
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
                return m_CodomainVar.ToList().AsReadOnly();
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
                return m_ParameterVar.ToList().AsReadOnly();
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


        ITemporalOperator m_TemporalOperator;

        /// <summary>
        /// %
        /// </summary>
        public ITemporalOperator TemporalOperator {
            get {
                return m_TemporalOperator;
            }
            set {
                if (IsCommitted)
                    throw new NotSupportedException("Not allowed to change after operator is committed.");
                m_TemporalOperator = value;
            }
        }

        MyDict m_FreeMeanValue;

        /// <summary>
        /// Notifies the solver that the mean value for a specific value is floating.
        /// An example is e.g. the pressure in the incompressible Navier-Stokes equation with all-walls boundary condition.
        /// - key: the name of some domain variable
        /// - value: false, if the mean value of the solution  is defined, true if the mean value  of the solution is floating (i.e. for some solution u, u + constant is also a solution).
        /// </summary>
        public IDictionary<string, bool> FreeMeanValue {
            get {
                if (m_FreeMeanValue == null) {
                    m_FreeMeanValue = new MyDict(this);
                }
                return m_FreeMeanValue;
            }
        }


        /// <summary>
        /// I hate shit like this class - so many dumb lines of code.
        /// </summary>
        class MyDict : IDictionary<string, bool> {

            XDifferentialOperatorMk2 owner;

            public MyDict(XDifferentialOperatorMk2 __owner) {
                owner = __owner;
                InternalRep = new Dictionary<string, bool>();
                foreach (string domName in __owner.DomainVar) {
                    InternalRep.Add(domName, false);
                }
            }

            Dictionary<string, bool> InternalRep;


            public bool this[string key] {
                get {
                    return InternalRep[key];
                }
                set {
                    if (!InternalRep.ContainsKey(key))
                        throw new ArgumentException("Must be a name of some domain variable.");
                    if (owner.IsCommitted)
                        throw new NotSupportedException("Changing is not allowed after operator is committed.");
                    InternalRep[key] = value;
                }
            }

            public ICollection<string> Keys => InternalRep.Keys;

            public ICollection<bool> Values => InternalRep.Values;

            public int Count => InternalRep.Count;

            public bool IsReadOnly => owner.IsCommitted;


            public void Add(string key, bool value) {
                throw new NotSupportedException("Addition/Removal of keys is not supported.");
            }

            public void Add(KeyValuePair<string, bool> item) {
                throw new NotSupportedException("Addition/Removal of keys is not supported.");
            }

            public void Clear() {
                throw new NotSupportedException("Addition/Removal of keys is not supported.");
            }

            public bool Contains(KeyValuePair<string, bool> item) {
                return InternalRep.Contains(item);
            }

            public bool ContainsKey(string key) {
                return InternalRep.ContainsKey(key);
            }


            public void CopyTo(KeyValuePair<string, bool>[] array, int arrayIndex) {
                (InternalRep as ICollection<KeyValuePair<string, bool>>).CopyTo(array, arrayIndex);
            }


            public IEnumerator<KeyValuePair<string, bool>> GetEnumerator() {
                throw new NotImplementedException();
            }

            public bool Remove(string key) {
                throw new NotSupportedException("Addition/Removal of keys is not supported.");
            }

            public bool Remove(KeyValuePair<string, bool> item) {
                throw new NotSupportedException("Addition/Removal of keys is not supported.");
            }

            public bool TryGetValue(string key, out bool value) {
                return InternalRep.TryGetValue(key, out value);
            }

            IEnumerator IEnumerable.GetEnumerator() {
                return InternalRep.GetEnumerator();
            }
        }

        #endregion


    }

}
