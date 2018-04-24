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
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Quadrature.FluxQuadCommon;
using BoSSS.Foundation.Quadrature.Linear;
using BoSSS.Platform;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Classic;
using System.Diagnostics;

namespace BoSSS.Foundation {

    /// <summary>
    /// This class represents a spatial operator which maps
    /// from (DG-) variables in the domain, identified and ordered by <see cref="DomainVar"/>
    /// to variables in the co-domain, identified and ordered by <see cref="CodomainVar"/>.
    /// </summary>
    public class SpatialOperator {

        /// <summary>
        /// Options for the treatment of edges at the boundary of a
        /// <see cref="SubGrid"/>.
        /// </summary>
        public enum SubGridBoundaryModes {

            /// <summary>
            /// Treats the edge as a real boundary edge (i.e., will use the
            /// boundary conditions defined by the flux function)
            /// </summary>
            BoundaryEdge = 0,

            /// <summary>
            /// Treats the edge as a standard inner edge (i.e., the inner edge
            /// flux will be evaluated as for any other inenr edge)
            /// </summary>
            InnerEdge,

            /// <summary>
            /// Treats the edge as an <i>open boundary</i> of the domain (i.e.,
            /// the inner edge flux will be used, but using the same values for
            /// the variables at both sides of the edge).
            /// </summary>
            OpenBoundary,

            /// <summary>
            /// Treats the edge as a standard inner edge (i.e., the inner edge
            /// flux will be evaluated as for any other inenr edge), but saves 
            /// also the fluxes across these edges (needed for LocalTimeStepping)
            /// </summary>
            InnerEdgeLTS,
        }


        /// <summary>
        /// Function Mapping from Domain Variable Degrees, Parameter Degrees and CoDomain Variable Degrees to the Quadrature Order
        /// </summary>
        public Func<int[], int[], int[], int> QuadOrderFunction;

        static string[] GetSubarray(string[] A, int i0, int len) {
            string[] r = new string[len];
            Array.Copy(A, i0, r, 0, len);
            return r;
        }

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="NoOfDomFields">number of domain fields</param>
        /// <param name="NoOfCodomFields">number of codomain fields</param>
        /// <param name="QuadOrderFunc">E.g., one of the members of <see cref="QuadOrderFunc"/>.</param>
        /// <param name="__varnames">
        /// names of domain variables 
        /// (entries 0 to (<paramref name="NoOfDomFields"/>-1)),
        /// followed by the names of the codomain variables
        /// (entries <paramref name="NoOfDomFields"/> to (<paramref name="NoOfDomFields"/>+<paramref name="NoOfCodomFields"/>-1));
        /// </param>
        public SpatialOperator(int NoOfDomFields, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, params string[] __varnames)
            : this(GetSubarray(__varnames, 0, NoOfDomFields), GetSubarray(__varnames, NoOfDomFields, NoOfCodomFields), QuadOrderFunc) {
            if(NoOfCodomFields + NoOfDomFields != __varnames.Length)
                throw new ArgumentException("mismatch between number of provided variable names and given number of domain and codomain fields.");
        }

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="NoOfDomFields"></param>
        /// <param name="NoOfCodomFields"></param>
        /// <param name="NoOfParameters"></param>
        /// <param name="QuadOrderFunc">E.g., one of the members of <see cref="QuadOrderFunc"/>.</param>
        /// <param name="__varnames">
        /// names of domain variables, followed by the names of the parameter variables,
        /// followed by the names of the codomain variables;
        /// </param>
        public SpatialOperator(int NoOfDomFields, int NoOfParameters, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, params string[] __varnames)
            : this(GetSubarray(__varnames, 0, NoOfDomFields), GetSubarray(__varnames, NoOfDomFields, NoOfParameters), GetSubarray(__varnames, NoOfDomFields + NoOfParameters, NoOfCodomFields), QuadOrderFunc) {
            if(NoOfCodomFields + NoOfDomFields + NoOfParameters != __varnames.Length)
                throw new ArgumentException("mismatch between number of provided variable names and given number of domain, parameter and codomain fields.");

        }

        /// <summary>
        /// constructor; 
        /// </summary>
        /// <param name="__DomainVar">
        /// variable names in the Domain of the spatial differential operator
        /// </param>
        /// <param name="__CoDomainVar">
        /// variable names in the Codomain of the spatial differential operator
        /// </param>
        /// <param name="QuadOrderFunc">E.g., one of the members of <see cref="QuadOrderFunc"/>.</param>
        public SpatialOperator(IList<string> __DomainVar, IList<string> __CoDomainVar, Func<int[], int[], int[], int> QuadOrderFunc)
            : this(__DomainVar, null, __CoDomainVar, QuadOrderFunc) {
        }

        /// <summary>
        /// constructor; 
        /// </summary>
        /// <param name="__DomainVar">
        /// variable names in the Domain of the spatial differential operator
        /// </param>
        /// <param name="__ParameterVar">
        /// variable names of the parameter variables; this list is optional, i.e. it 
        /// can be null;
        /// </param>
        /// <param name="__CoDomainVar">
        /// variable names in the Codomain of the spatial differential operator
        /// </param>
        /// <param name="QuadOrderFunc">E.g., one of the members of <see cref="QuadOrderFunc"/>.</param>
        public SpatialOperator(IList<string> __DomainVar, IList<string> __ParameterVar, IList<string> __CoDomainVar, Func<int[], int[], int[], int> QuadOrderFunc) {
            m_DomainVar = new string[__DomainVar.Count];
            for(int i = 0; i < m_DomainVar.Length; i++) {
                if(Array.IndexOf<string>(m_DomainVar, __DomainVar[i]) >= 0)
                    throw new ArgumentException("error in domain variables list; identifier \"" + __DomainVar[i] + "\" appears twice.", "__DomainVar");
                m_DomainVar[i] = __DomainVar[i];
            }

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

            m_EquationComonents = new SortedList<string, List<IEquationComponent>>(__CoDomainVar.Count);
            foreach(var f in __CoDomainVar) {
                m_EquationComonents.Add(f, new List<IEquationComponent>());
            }
            m_EquationComponentsHelper = new _EquationComponents(this);
            this.QuadOrderFunction = QuadOrderFunc;
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
        private void Verify() {
            foreach(var comps in m_EquationComonents.Values) {
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
        /// Evaluates this operator for given DG fields;
        /// </summary>
        /// <param name="DomainMapping">
        /// the domain variables, or "input data" for the operator; the number
        /// of elements must be equal to the number of elements in the
        /// <see cref="DomainVar"/>-list;
        /// </param>
        /// <param name="Params">
        /// List of parameter fields; May be null
        /// </param>
        /// <param name="CodomainMapping">
        /// the co-domain variables, or "output" for the evaluation of the
        /// operator; the number of elements must be equal to the number of
        /// elements in the <see cref="CodomainVar"/>-list;
        /// </param>
        /// <param name="alpha">
        /// scaling of the operator 
        /// </param>
        /// <param name="beta">
        /// scaling of the accumulator (<paramref name="CodomainMapping"/>);
        /// </param>
        /// <param name="sgrd">
        /// subgrid, for restricted evaluation; null indicates evaluation on
        /// the full grid.
        /// </param>
        /// <param name="bndMode">
        /// Treatment of subgrid boundaries, if <paramref name="sgrd"/> is not
        /// null. See <see cref="Evaluator"/>
        /// </param>
        /// <param name="qInsEdge">
        /// Optional definition of the edge quadrature scheme. Since this
        /// already implies a domain of integration, must be null if
        /// <paramref name="sgrd"/> is not null.
        /// </param>
        /// <param name="qInsVol">
        /// Optional definition of the volume quadrature scheme. Since this
        /// already implies a domain of integration, must be null if
        /// <paramref name="sgrd"/> is not null.
        /// </param>
        /// <remarks>
        /// If some of the input data, <paramref name="DomainMapping"/>, is
        /// contained in the  output data, <paramref name="CodomainMapping"/>,
        /// these DG fields will be cloned to ensure correct operation of the
        /// operator evaluation.<br/>
        /// It is not a good choice to use this function if this operator
        /// should be evaluated multiple times and  contains linear components
        /// (i.e. <see cref="ContainsLinear"/> returns true); If the latter is
        /// the case, the matrix which represents the linear components of the
        /// operator must be computed first, which is computational- and
        /// memory-intensive; After execution of this method, the matrix will
        /// be lost; If multiple evaluation is desired, the
        /// <see cref="Evaluator"/>-class should be used, in which the matrix
        /// of the operator will persist; However, if no linear components are
        /// present, the performance of this function should be almost
        /// comparable to the use of the <see cref="Evaluator"/>-class;
        /// </remarks>
        public void Evaluate(double alpha, double beta,
                             CoordinateMapping DomainMapping, IList<DGField> Params, CoordinateMapping CodomainMapping,
                             SubGrid sgrd = null,
                             EdgeQuadratureScheme qInsEdge = null,
                             CellQuadratureScheme qInsVol = null,
                             SpatialOperator.SubGridBoundaryModes bndMode = SubGridBoundaryModes.OpenBoundary,
                             double time = double.NaN) //
        {
            using(new FuncTrace()) {
                if(sgrd != null && (qInsEdge != null || qInsVol != null))
                    throw new ArgumentException("Specification of Subgrid and quadrature schemes is exclusive: not allowed to specify both at the same time.", "sgrd");
#if DEBUG
                Verify(); //this has already been done during this.Commit()
#endif

                IList<DGField> _DomainFields = DomainMapping.Fields;
                IList<DGField> _CodomainFields = CodomainMapping.Fields;
                DGField[] _DomainFieldsRevisited = new DGField[_DomainFields.Count];

                bool a = false;
                for(int i = 0; i < _DomainFields.Count; i++) {
                    DGField f = _DomainFields[i];

                    if(_CodomainFields.Contains(f)) {
                        // some of the domain variables (input data)
                        // is also a member of the codomain variables (output);
                        // the data need to be cloned to provide correct results
                        a = true;
                        _DomainFieldsRevisited[i] = (DGField)f.Clone();
                    } else
                        _DomainFieldsRevisited[i] = f;
                }

                CoordinateMapping domainMappingRevisited;
                if(a)
                    domainMappingRevisited = new CoordinateMapping(_DomainFieldsRevisited);
                else
                    domainMappingRevisited = DomainMapping;

                if(sgrd != null) {
                    CellMask cm = (sgrd == null) ? null : sgrd.VolumeMask;
                    EdgeMask em = (sgrd == null) ? null : sgrd.AllEdgesMask;

                    qInsEdge = new EdgeQuadratureScheme(true, em);
                    qInsVol = new CellQuadratureScheme(true, cm);
                }

                var ev = GetEvaluatorEx(
                    domainMappingRevisited, Params, CodomainMapping,
                    qInsEdge,
                    qInsVol);

                if(sgrd != null)
                    ev.ActivateSubgridBoundary(sgrd.VolumeMask, bndMode);
                CoordinateVector outp = new CoordinateVector(CodomainMapping);
                ev.time = time;
                ev.Evaluate<CoordinateVector>(alpha, beta, outp);
            }
        }

        /// <summary>
        /// Another wrapper for
        /// <see cref="Evaluate(double,double,CoordinateMapping,IList{DGField},CoordinateMapping,SubGrid,EdgeQuadratureScheme,CellQuadratureScheme,SpatialOperator.SubGridBoundaryModes)"/>
        /// </summary>
        /// <param name="DomainFields">
        /// the domain variables; the number of elements
        /// must be equal to the number of elements in the <see cref="DomainVar"/>-list;
        /// </param>
        /// <param name="CodomainFields">
        /// the codomain variables; the number of elements
        /// must be equal to the number of elements in the <see cref="CodomainVar"/>-list;
        /// </param>
        /// <remarks>
        /// If some of the input field, i.e. some element of <paramref name="DomainFields"/>, is contained in the 
        /// output fields, i.e. in the list <paramref name="CodomainFields"/>, these DG fields will be cloned to ensure 
        /// correct operation of the operator evaluation.<br/>
        /// It is not a good choice to use this 
        /// function if this operator should be evaluated multiple times
        /// and 
        /// contains linear components (i.e. <see cref="ContainsLinear"/> returns true);
        /// If the later one is the case, the matrix which represents
        /// the linear components of the operator must be computed first, which is computational- and
        /// memory - intensive;
        /// After execution of this method, the matrix will be lost;
        /// If multiple evaluation is desired, the <see cref="Evaluator"/>-class should be used,
        /// in which the matrix of the operator will persist;
        /// However, if no linear components are present, the performance of this function should
        /// be almost comparable to the use of the <see cref="Evaluator"/>-class;
        /// </remarks>
        public void Evaluate(IList<DGField> DomainFields, IList<DGField> CodomainFields, double time = double.NaN) {
            if(DomainFields.Count != m_DomainVar.Length)
                throw new ArgumentException("wrong number of domain fields", "DomainFields");
            if(CodomainFields.Count != m_CodomainVar.Length)
                throw new ArgumentException("wrong number of domain fields", "CodomainFields");

            CoordinateMapping inp = new CoordinateMapping(DomainFields);
            CoordinateMapping outp = new CoordinateMapping(CodomainFields);

            this.Evaluate(1.0, 0.0, inp, null, outp, null, null, null, SubGridBoundaryModes.OpenBoundary, time);
        }

        /// <summary>
        /// And another wrapper.
        /// </summary>
        public void Evaluate(params DGField[] f) {
            this.Evaluate(double.NaN, f);
        }


        /// <summary>
        /// And another wrapper.
        /// </summary>
        public void Evaluate(double time, params DGField[] f) {
            this.Evaluate(time, null, SubGridBoundaryModes.OpenBoundary, f);
        }

        /// <summary>
        /// And another wrapper.
        /// </summary>
        public void Evaluate(double time, SubGrid subGrid = null, SubGridBoundaryModes subGridBoundaryMode = SubGridBoundaryModes.OpenBoundary, params DGField[] f) {
            if(DomainVar.Count + ParameterVar.Count + CodomainVar.Count != f.Length)
                throw new ArgumentException("wrong number of domain/parameter/codomain fields", "f");



            CoordinateMapping inp = new CoordinateMapping(f.GetSubVector(0, this.DomainVar.Count));
            DGField[] Parameters = f.GetSubVector(this.DomainVar.Count, this.ParameterVar.Count);
            CoordinateMapping outp = new CoordinateMapping(f.GetSubVector(this.DomainVar.Count + this.ParameterVar.Count, this.CodomainVar.Count));

            this.Evaluate(1.0, 0.0, inp, Parameters, outp, subGrid, null, null, subGridBoundaryMode, time);
        }


        /// <summary>
        /// another wrapper for
        /// <see cref="Evaluate(double,double,CoordinateMapping,IList{DGField},CoordinateMapping,SubGrid,EdgeQuadratureScheme,CellQuadratureScheme,SpatialOperator.SubGridBoundaryModes)"/>
        /// </summary>
        /// <param name="NoOfDomainVar"></param>
        /// <param name="NoOfCodomainVar"></param>
        /// <param name="fields">
        /// a list of <paramref name="NoOfDomainVar"/>+<paramref name="NoOfCodomainVar"/> fields;
        /// </param>
        public void Evaluate(uint NoOfDomainVar, uint NoOfCodomainVar, double time, params DGField[] fields) {
            if(fields.Length != (NoOfDomainVar + NoOfCodomainVar))
                throw new ArgumentException("wrong number of fields.", "fields");

            DGField[] DomF = new DGField[NoOfDomainVar];
            Array.Copy(fields, 0, DomF, 0, NoOfDomainVar);
            DGField[] CoDomF = new DGField[NoOfCodomainVar];
            Array.Copy(fields, NoOfDomainVar, CoDomF, 0, NoOfCodomainVar);

            Evaluate(DomF, CoDomF, time);
        }


        /// <summary>
        /// simplified version of 
        /// <see cref="ComputeMatrixEx{M,V}"/>;
        /// </summary>
        public void ComputeMatrix<M, V>(
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
            M Matrix, V AffineOffset, bool OnlyAffine, double time = 0.0,
            SubGrid sgrd = null)
            where M : IMutableMatrixEx
            where V : IList<double> {

            var GridDat = CheckArguments(DomainMap, Parameters, CodomainMap);


            var ev = this.GetMatrixBuilder(DomainMap, Parameters, CodomainMap,
                new EdgeQuadratureScheme(true, sgrd == null ? null : sgrd.AllEdgesMask),
                 new CellQuadratureScheme(true, sgrd == null ? null : sgrd.VolumeMask));

            ev.time = time;

            if(OnlyAffine) {
                ev.ComputeAffine(AffineOffset);
            } else {
                ev.ComputeMatrix(Matrix, AffineOffset);
            }
        }

        public int GetOrderFromQuadOrderFunction(UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap) {
            /// Compute Quadrature Order
            int order;
            int[] DomainDegrees = DomainMap.BasisS.Select(f => f.Degree).ToArray();
            int[] CodomainDegrees = CodomainMap.BasisS.Select(f => f.Degree).ToArray();
            int[] ParameterDegrees;
            if(Parameters != null && Parameters.Count != 0) {
                ParameterDegrees = Parameters.Select(f => f == null ? 0 : f.Basis.Degree).ToArray();
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
        /// simplified version of 
        /// <see cref="ComputeMatrixEx{M,V}"/>;
        /// </summary>
        public MsrMatrix ComputeMatrix(
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap, double time = 0.0, SubGrid SubGrid = null) {


            int RowBlkSize = (CodomainMap.MaxTotalNoOfCoordinatesPerCell == CodomainMap.MinTotalNoOfCoordinatesPerCell) ? CodomainMap.MaxTotalNoOfCoordinatesPerCell : 1;
            int ColBlkSize = (DomainMap.MaxTotalNoOfCoordinatesPerCell == DomainMap.MinTotalNoOfCoordinatesPerCell) ? DomainMap.MaxTotalNoOfCoordinatesPerCell : 1;

            MsrMatrix Matrix = new MsrMatrix(CodomainMap.LocalLength, (int)DomainMap.GlobalCount, RowBlkSize, ColBlkSize);

            double[] dummyRHS = new double[Matrix.RowPartitioning.LocalLength];
            ComputeMatrix(
                DomainMap, Parameters, CodomainMap,
                Matrix, dummyRHS,
                false, time, SubGrid);

            return Matrix;
        }

        /*
        /// <summary>
        /// Simplified version of <see cref="ComputeAffine{V}"/>.
        /// </summary>
        virtual public double[] ComputeAffine(UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap, double time = 0.0) {
            double[] affine = new double[CodomainMap.LocalLength];

            this.ComputeAffine(DomainMap, Parameters, CodomainMap, affine, OnlyBoundaryEdges: false, time: time);

            return affine;
        }
        */

        /*
        /// <summary>
        /// computes the affine offset and/or matrix of the operator, expert
        /// version;
        /// </summary>
        /// <param name="DomainMap">
        /// the mapping which is used to compute column indices into
        /// <paramref name="Matrix"/>;
        /// </param>
        /// <param name="Parameters">
        /// The parameter variables (of this differential operator);
        /// The number of elements in the list must match the parameter count
        /// of the differential operator (see
        /// <see cref="SpatialOperator.ParameterVar"/>);  It is allowed to set
        /// an entry to 'null', in this case the values of the parameter field
        /// are assumed to be 0.0; If the differential operator contains no
        /// parameters, this argument can be null;
        /// </param>
        /// <param name="CodomainMap">
        /// the mapping which is used to compute row indices into
        /// <paramref name="Matrix"/> and <paramref name="AffineOffset"/>.
        /// </param>
        /// <param name="Matrix">
        /// Acc output: the matrix which represents the linear part of this
        /// operator, according to the mapping given by
        /// <paramref name="DomainMap"/> and <paramref name="CodomainMap"/>,
        /// is <b>ACCUMULATED</b> here; <br/>
        /// Setting all matrix entries to 0.0 is left to the user;
        /// </param>
        /// <param name="AffineOffset">
        /// Acc output: the vector which represents the affine part of this
        /// operator, according to the mapping given by
        /// <paramref name="DomainMap"/> and <paramref name="CodomainMap"/>,
        /// is <b>ACCUMULATED</b> here; <br/>
        /// Setting all vector entries to 0.0 is left to the user;
        /// </param>
        /// <param name="OnlyAffine">
        /// If true, only the <paramref name="AffineOffset"/> is computed, and
        /// the <paramref name="Matrix"/> is not touched (can be null);
        /// </param>
        /// <remarks>
        /// The operator assembly must be finalized before by calling
        /// <see cref="Commit"/> before this method can be called.
        /// </remarks>
        /// <param name="edgeRule">
        /// Quadrature rule and domain for edge integration; specifying this is exclusive with <paramref name="edgeQuadScheme"/>, i.e. both cannot be unequal null at the same time.
        /// </param>
        /// <param name="edgeQuadScheme">
        /// Quadrature scheme for edge integration; specifying this is exclusive with <paramref name="edgeRule"/>, i.e. both cannot be unequal null at the same time.
        /// </param>
        /// <param name="volRule">
        /// Quadrature rule and domain for volume integration; specifying this is exclusive with <paramref name="volQuadScheme"/>, i.e. both cannot be unequal null at the same time.
        /// </param>
        /// <param name="volQuadScheme">
        /// Quadrature scheme for volume integration; specifying this is exclusive with <paramref name="volRule"/>, i.e. both cannot be unequal null at the same time.
        /// </param>
        /// <param name="SubGridBoundaryMask">
        /// </param>
        /// <param name="ParameterMPIExchange">
        /// Determines whether parameter fields have to exchange ghost cell
        /// data before the assembly of the operator.
        /// </param>
        /// <param name="time"></param>
        virtual public void ComputeMatrixEx<M, V>(
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
            M Matrix, V AffineOffset,bool OnlyAffine = false,
            double time = 0.0, 
            EdgeQuadratureScheme edgeQuadScheme = null, CellQuadratureScheme volQuadScheme = null,
            ICompositeQuadRule<QuadRule> edgeRule = null, ICompositeQuadRule<QuadRule> volRule = null,
            BitArray SubGridBoundaryMask = null,
            bool ParameterMPIExchange = true)
            where M : IMutableMatrix
            where V : IList<double> //
        {
            //

            if(edgeQuadScheme != null && edgeRule != null)
                throw new ArgumentException();
            if(volQuadScheme != null && volRule != null)
                throw new ArgumentException();

            //if (edgeQrCtx == null)
            //    edgeQrCtx = new EdgeQuadratureScheme(true);

            var GridDat = CheckArguments(DomainMap, Parameters, CodomainMap);

            int order = 0;
            if(edgeRule == null || volRule == null) 
                order = this.GetOrderFromQuadOrderFunction(DomainMap, Parameters, CodomainMap);
            
            if(edgeRule == null)
                edgeRule = edgeQuadScheme.SaveCompile(GridDat, order);
            if(volRule == null)
                volRule = volQuadScheme.SaveCompile(GridDat, order);

            Internal_ComputeMatrixEx(GridDat, DomainMap, Parameters, CodomainMap, Matrix, AffineOffset, OnlyAffine, time,
                edgeRule,
                volRule,
                SubGridBoundaryMask, ParameterMPIExchange);
        }
        */

        /*
        /// <summary>
        /// An important special case of <see cref="ComputeMatrixEx{M,V}"/>, in order to compute only
        /// the affine offset.
        /// </summary>
        /// <typeparam name="V"></typeparam>
        /// <param name="DomainMap">Domain mapping - needs to be given to determine quadrature order.</param>
        /// <param name="Parameters">Parameter variables</param>        
        /// <param name="CodomainMap">Codomain Basis</param>
        /// <param name="AffineOffset">Output</param>
        /// <param name="OnlyBoundaryEdges">
        /// If true, integration is only carried out on the boundary edges of the 
        /// computational domain. regarding the affine offset, this covers most application scenarios
        /// </param>
        /// <param name="edgeQr">
        /// Quadrature instruction for edge integrals: if specified (not null), <paramref name="OnlyBoundaryEdges"/> must be false.
        /// </param>
        /// <param name="volQr">
        /// Quadrature instruction for volume integrals: if specified (not null), <paramref name="OnlyBoundaryEdges"/> must be false.
        /// </param>
        /// <param name="time"></param>
        virtual public void ComputeAffine<V>(UnsetteledCoordinateMapping DomainMap,
            IList<DGField> Parameters,
            UnsetteledCoordinateMapping CodomainMap,
            V AffineOffset,
            bool OnlyBoundaryEdges = true, double time = 0.0,
            EdgeQuadratureScheme edgeQr = null, CellQuadratureScheme volQr = null)
            where V : IList<double> {

            var GridDat = CodomainMap.GridDat;
            if (Parameters != null)
                foreach (var prm in Parameters)
                    if (!object.ReferenceEquals(prm.GridDat, GridDat))
                        throw new ArgumentException(string.Format("parameter field {0} is assigned to a different grid.", prm.Identification));

            //Using order zero for DomainMap will lead to inconsistent (and possibly insufficient) quadrature order!!! 
            //UnsetteledCoordinateMapping DomainMap;
            //Basis b = new Basis(GridDat, 0);
            //Basis[] B = new Basis[this.DomainVar.Count];
            //B.SetAll(b);
            //DomainMap = new UnsetteledCoordinateMapping(B);


            if (OnlyBoundaryEdges) {
                if (edgeQr != null)
                    throw new ArgumentException("If 'OnlyBoundaryEdges == true', 'edgeQr' must be null!", "edgeQr");
                if (volQr != null)
                    throw new ArgumentException("If 'OnlyBoundaryEdges == true', 'volQr' must be null!", "volQr");

                volQr = new CellQuadratureScheme(true, CellMask.GetEmptyMask(GridDat));
                edgeQr = new EdgeQuadratureScheme(true, GridDat.GetBoundaryEdgeMask());
            }

            this.ComputeMatrixEx(
                DomainMap, Parameters, CodomainMap,
                default(MsrMatrix), AffineOffset,
                OnlyAffine: true, time:time,
                volQuadScheme: volQr, edgeQuadScheme: edgeQr);
        }
        */




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

            var comps = m_EquationComonents[CodomVar];
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
        SortedList<string, List<IEquationComponent>> m_EquationComonents;

        _EquationComponents m_EquationComponentsHelper;

        /// <summary>
        /// for each variable in <see cref="CodomainVar"/>, a
        /// collection of equation components that define the operator.
        /// 
        /// </summary>
        public _EquationComponents EquationComponents {
            get {
                return m_EquationComponentsHelper;
            }
        }

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
                return this.m_EquationComonents.Values.Sum(x => x.Count);
            }
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

        }


        /// <summary>
        /// implementation of <see cref="EquationComponents" />;
        /// </summary>
        public class _EquationComponents : IEnumerable<KeyValuePair<string, IEnumerable<IEquationComponent>>> {

            internal _EquationComponents(SpatialOperator owner) {
                m_owner = owner;
            }

            SpatialOperator m_owner;

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
                        return m_owner.m_EquationComonents[EqnName].AsReadOnly();
                    else
                        return m_owner.m_EquationComonents[EqnName];
                }
            }

            #region IEnumerable<KeyValuePair<string,IEnumerable<IEquationComponent>> Members

            /// <summary>
            /// Gets the enumerator over the equation components of the owner
            /// </summary>
            /// <returns>An enumerator</returns>
            public IEnumerator<KeyValuePair<string, IEnumerable<IEquationComponent>>> GetEnumerator() {
                return m_owner.m_EquationComonents.Select(
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
        /// returns a collection of equation components of a certain type (<typeparamref name="T"/>)
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="CatParams">
        /// if true, parameter variables (see <see cref="IEquationComponent.ParameterOrdering"/>)
        /// are concatenated with domain variable names (see <see cref="IEquationComponent."ArgumentOrdering/>).
        /// </param>
        /// <param name="F">
        /// optional filter;
        /// should return true, if the component should be added, false if not; 
        /// </param>
        /// <param name="vectorizer">
        /// vectorizer option: translate some equation component to another one
        /// </param>
        public EquationComponentArgMapping<T>[] GetArgMapping<T>(bool CatParams = false, Func<T, bool> F = null, Func<IEquationComponent, IEquationComponent> vectorizer = null) where T : IEquationComponent {
            if(!IsCommited)
                throw new ApplicationException("Commit() has to be called prior to this method.");

            int Gamma = CodomainVar.Count;

            var ret = new EquationComponentArgMapping<T>[Gamma];
            for(int i = 0; i < Gamma; i++) {
                var codName = this.m_CodomainVar[i];
                ret[i] = new EquationComponentArgMapping<T>(this,
                    codName,
                    this.m_DomainVar,
                    CatParams ? this.ParameterVar : null,
                    F, vectorizer);
            }

            return ret;
        }



        /// <summary>
        /// returns true, if this spatial differential operator contains any 
        /// linear component (i.e. an object that implements
        /// <see cref="ILinear2ndDerivativeFlux"/> or <see cref="ILinearSource"/> or <see cref="ILinearFlux"/>),
        /// otherwise returns false;
        /// </summary>
        /// <returns></returns>
        public bool ContainsLinear() {
            foreach(string var in m_CodomainVar) {
                if(ContainsLinearPerCodomainVar(var))
                    return true;
            }

            return false;
        }

        /// <summary>
        /// returns true, if this spatial differential operator contains any 
        /// nonlinear component (i.e. an object that implements
        /// <see cref="INonlinearFlux"/> or <see cref="INonlinearFluxEx"/> or <see cref="INonlinearSource"/>
        /// or <see cref="IDualValueFlux"/>),
        /// otherwise returns false;
        /// </summary>
        public bool ContainsNonlinear {
            get {
                foreach(string var in m_CodomainVar) {
                    if(ContainsNonlinearPerCodomainVar(var))
                        return true;
                }

                return false;
            }
        }

        /// <summary>
        /// True, if this operator contains any edge term.
        /// </summary>
        public bool RequiresEdgeQuadrature {
            get {
                foreach(string var in m_CodomainVar) {
                    if(ReqNonlinEdgeQuad_PerCodomainVar(var))
                        return true;
                }

                return false;
            }
        }

        /// <summary>
        /// True, if this operator contains any volume term.
        /// </summary>
        public bool RequiresVolumeQuadrature {
            get {
                foreach(string var in m_CodomainVar) {
                    if(ReqNonlVolumeQuad_PerCodomainVar(var))
                        return true;
                }

                return false;
            }
        }



        /// <summary>
        /// returns true, if any of the equation components associated with 
        /// variable <paramref name="CodomVar"/> is linear
        /// (i.e. an object that implements
        /// <see cref="ILinear2ndDerivativeFlux"/> or <see cref="ILinearSource"/> or <see cref="ILinearFlux"/>);
        /// These components are <see cref="EquationComponents"/>[<paramref name="CodomVar"/>];
        /// </summary>
        /// <param name="CodomVar">
        /// identifies the codomain variable name, must be a member of
        /// <see cref="CodomainVar"/>;
        /// </param>
        /// <returns></returns>
        public bool ContainsLinearPerCodomainVar(string CodomVar) {
            if(Array.IndexOf<string>(m_CodomainVar, CodomVar) < 0)
                throw new ArgumentException("given codomain variable name is not known in this spatial differential operator", "CodomVar");

            foreach(object o in EquationComponents[CodomVar]) {
                Type[] interfaces = o.GetType().GetInterfaces();

                //if (Array.IndexOf<Type>(interfaces, typeof(ILinear2ndDerivativeFlux)) >= 0)
                //    return true;
                //if (Array.IndexOf<Type>(interfaces, typeof(ILinearDualValueFlux)) >= 0)
                //    return true;
                //if (Array.IndexOf<Type>(interfaces, typeof(ILinearSource)) >= 0)
                //    return true;
                //if (Array.IndexOf<Type>(interfaces, typeof(ILinearDerivativeSource)) >= 0)
                //    return true;
                //if (Array.IndexOf<Type>(interfaces, typeof(ILinearFlux)) >= 0)
                //    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IVolumeForm_GradUxGradV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeform_GradUxV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IVolumeForm_GradUxGradV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeform_UxGradV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeSource_V)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeSource_GradV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeform_UxV)) >= 0)
                    return true;
            }

            return false;
        }

        /// <summary>
        /// returns true, if any of the equation components associated with 
        /// variable <paramref name="CodomVar"/> is nonlinear
        /// </summary>
        /// <param name="CodomVar">
        /// identifies the codomain variable name, must be a member of
        /// <see cref="CodomainVar"/>;
        /// </param>
        /// <returns></returns>
        public bool ContainsNonlinearPerCodomainVar(string CodomVar) {
            return (ReqNonlinEdgeQuad_PerCodomainVar(CodomVar) || ReqNonlVolumeQuad_PerCodomainVar(CodomVar));
        }

        /// <summary>
        /// returns true, if any of the equation components associated with 
        /// variable <paramref name="CodomVar"/> contains nonlinear integrands
        /// on edges (these are objects that implement <see cref="INonlinearFlux"/>,
        /// <see cref="INonlinearFluxEx"/> or <see cref="INonlinearSource"/>);
        /// </summary>
        /// <param name="CodomVar">
        /// identifies the codomain variable name, must be a member of
        /// <see cref="CodomainVar"/>;
        /// </param>
        /// <returns></returns>
        public bool ReqNonlVolumeQuad_PerCodomainVar(string CodomVar) {
            if(Array.IndexOf<string>(m_CodomainVar, CodomVar) < 0)
                throw new ArgumentException("given codomain variable name is not known in this spatial differential operator", "CodomVar");

            foreach(object o in EquationComponents[CodomVar]) {
                Type[] interfaces = o.GetType().GetInterfaces();

                if(Array.IndexOf<Type>(interfaces, typeof(IVolumeForm)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(INonlinearFlux)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(INonlinearFluxEx)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(INonlinearSource)) >= 0)
                    return true;
            }

            return false;
        }

        /// <summary>
        /// returns true, if any of the equation components associated with 
        /// variable <paramref name="CodomVar"/> contains nonlinear integrands
        /// on edges (these are objects that implement <see cref="INonlinearFlux"/>,
        /// <see cref="INonlinearFluxEx"/> or <see cref="IDualValueFlux"/>);
        /// </summary>
        /// <param name="CodomVar">
        /// identifies the codomain variable name, must be a member of
        /// <see cref="CodomainVar"/>;
        /// </param>
        /// <returns></returns>
        bool ReqNonlinEdgeQuad_PerCodomainVar(string CodomVar) {
            if(Array.IndexOf<string>(m_CodomainVar, CodomVar) < 0)
                throw new ArgumentException("given codomain variable name is not known in this spatial differential operator", "CodomVar");

            foreach(object o in EquationComponents[CodomVar]) {
                Type[] interfaces = o.GetType().GetInterfaces();

                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeForm)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(INonlinearFlux)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(INonlinearFluxEx)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IDualValueFlux)) >= 0)
                    return true;
            }

            return false;
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


        /// <summary>
        /// constructs a new Instance of <see cref="Evaluator"/>
        /// </summary>
        /// <returns></returns>
        /// <remarks>
        /// The
        /// operator assembly must be finalized before by calling <see cref="Commit"/> before this method can be called.
        /// </remarks>
        /// <param name="CodomainVarMap">
        /// used to compute indices into the result vector
        /// </param>
        /// <param name="DomainFields">
        /// domain which are evaluated to compute fluxes, ...
        /// </param>
        /// <param name="ParameterMap">
        /// The parameter variables (of this differential operator);
        /// The number of elements in the list must match the parameter count of the differential operator
        /// (see <see cref="SpatialOperator.ParameterVar"/>);
        /// It is allowed to set an entry to 'null', in this case the values of the parameter field
        /// are assumed to be 0.0;
        /// If the differential operator contains no parameters, this argument can be null;
        /// </param>
        /// <param name="edgeQrCtx">optional quadrature instruction for edges</param>
        /// <param name="volQrCtx">optional quadrature instruction for volumes/cells</param>
        /// <param name="subGridBoundaryTreatment">
        /// Optional definition of the treatment of edges at the boundary of a
        /// <see cref="SubGrid"/>. By default, they will be treated as boundary
        /// edges (see
        /// <see cref="SubGridBoundaryModes.BoundaryEdge"/>)
        /// </param>
        /// <param name="sgrd">
        /// </param>
        public virtual IEvaluatorNonLin GetEvaluatorEx(
            IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
            EdgeQuadratureScheme edgeQrCtx = null,
            CellQuadratureScheme volQrCtx = null) //
        {

            using(new FuncTrace()) {

                /// This is already done in the constructor of Evaluator
#if DEBUG
                if(!m_IsCommited)
                    throw new ApplicationException("operator assembly must be finalized before by calling 'Commit' before this method can be called.");
#endif


                var e = new EvaluatorNonLin(this, new CoordinateMapping(DomainFields), ParameterMap, CodomainVarMap, edgeQrCtx, volQrCtx);

                return e;
            }
        }

        public virtual IEvaluatorLinear GetMatrixBuilder(
            UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
            EdgeQuadratureScheme edgeQrCtx = null,
            CellQuadratureScheme volQrCtx = null) //
        {

            using(new FuncTrace()) {

                /// This is already done in the constructor of Evaluator
#if DEBUG
                if(!m_IsCommited)
                    throw new ApplicationException("operator assembly must be finalized before by calling 'Commit' before this method can be called.");
#endif


                var e = new EvaluatorLinear(this, DomainVarMap, ParameterMap, CodomainVarMap, edgeQrCtx, volQrCtx);

                return e;
            }
        }




        /// <summary>
        /// Container for the evaluation of nonlinear fluxes/sources
        /// </summary>
        abstract public class EvaluatorBase : IEvaluator {

            SpatialOperator m_Owner;

            /// <summary>
            /// the operator used to construct this object
            /// </summary>
            public SpatialOperator Owner {
                get {
                    return m_Owner;
                }
            }

            /// <summary>
            /// ctor
            /// </summary>
            protected internal EvaluatorBase(
                SpatialOperator owner,
                UnsetteledCoordinateMapping DomainVarMap,
                IList<DGField> ParameterMap,
                UnsetteledCoordinateMapping CodomainVarMap) //
            {
                using (var tr = new FuncTrace()) {
                    MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                    if (DomainVarMap.NoOfVariables != owner.DomainVar.Count) {
                        throw new ArgumentException("wrong number of domain variables provided.");
                    }
                    this.m_Parameters = new DGField[owner.ParameterVar.Count];
                    if (CodomainVarMap.NoOfVariables != owner.CodomainVar.Count) {
                        throw new ArgumentException("wrong number of codomain variables provided.");
                    }

                    IEnumerable<Basis> allBasis = DomainVarMap.BasisS;
                    if (ParameterMap != null) {
                        allBasis = allBasis.Union(ParameterMap.Select(f => f.Basis));
                    }
                    allBasis = allBasis.Union(CodomainVarMap.BasisS);
                    IGridData grdDat = allBasis.First().GridDat;
                    foreach (var b in allBasis) {
                        if (!object.ReferenceEquals(grdDat, b.GridDat)) {
                            throw new ArgumentException("all fields (domain, parameter, codomain) must be defined on the same grid.");
                        }
                    }

                    m_Owner = owner;
                    m_CodomainMapping = CodomainVarMap;
                    m_DomainMapping = DomainVarMap;
                    m_Parameters = (ParameterMap != null) ? ParameterMap.ToArray() : new DGField[0];

                    if (!m_Owner.IsCommited)
                        throw new ApplicationException("operator assembly must be finalized before by calling 'Commit' before this method can be called.");

                    order = owner.GetOrderFromQuadOrderFunction(m_DomainMapping, ParameterMap, CodomainVarMap);

                    m_OperatorCoefficients = new CoefficientSet() {
                        CellLengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Cells.cj,
                        EdgeLengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Edges.h_min_Edge,
                        UserDefinedValues = new Dictionary<string, object>()
                    };
                }
            }

            CoefficientSet m_OperatorCoefficients;

            /// <summary>
            /// Stuff passed to equation components which implement <see cref="IEquationComponentCoefficient"/>.
            /// </summary>
            virtual public CoefficientSet OperatorCoefficients  {
                get {
                    return m_OperatorCoefficients;
                }
                set {
                    m_OperatorCoefficients = value;
                }
            }


            /// <summary>
            /// Sets the coefficients for all equation components of the operator which implement <see cref="IEquationComponentCoefficient"/>.
            /// </summary>
            protected void UpdateCoefficients() {
                int[] DomDGdeg = this.DomainMapping.BasisS.Select(b => b.Degree).ToArray();
                int[] CodDGdeg = this.CodomainMapping.BasisS.Select(b => b.Degree).ToArray();
                string[] DomNames = m_Owner.DomainVar.ToArray();
                string[] CodNames = m_Owner.CodomainVar.ToArray();


                Debug.Assert(CodNames.Length == CodDGdeg.Length);
                for(int iCod = 0; iCod < CodDGdeg.Length; iCod++) {
                    var comps = m_Owner.m_EquationComonents[CodNames[iCod]];
                    foreach(var c in comps) {
                        if(c is IEquationComponentCoefficient) {
                            var ce = c as IEquationComponentCoefficient;

                            int[] DomDGdeg_cd = new int[ce.ArgumentOrdering.Count];
                            for(int i = 0; i < DomDGdeg_cd.Length; i++) {
                                string domName = ce.ArgumentOrdering[i];
                                int idx = Array.IndexOf(DomNames, domName);
                                DomDGdeg_cd[i] = DomDGdeg[idx];
                            }


                            ce.CoefficientUpdate(m_OperatorCoefficients, DomDGdeg_cd, CodDGdeg[iCod]);
                        }
                    }
                }
            }


            /// <summary>
            /// Quadrature order to compile quadrature schemes
            /// </summary>
            public int order {
                get;
                private set;
            }


            /*
            public void SetParameters(IEnumerable<DGField> pFields) {
                if (pFields.Count() != m_Owner.ParameterVar.Count) {
                    throw new ArgumentException("wrong number of parameter variables provided.");
                }

                foreach(var p in pFields) {
                    if(p!= null) {
                        if(!object.ReferenceEquals(p.GridDat, this.GridData))
                            throw new ArgumentException("Mismatch between parameter field grid and grid for the evaluator.");
                    }
                }

                m_Parameters = pFields.ToArray();
            }
            */



            


            SubGridBoundaryModes m_SubGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge;

            /// <summary>
            /// 
            /// </summary>
            protected CellMask m_SubGrid_InCells;


            public SubGridBoundaryModes SubGridBoundaryTreatment {
                get {
                    return m_SubGridBoundaryTreatment;
                }
            }


            public void ActivateSubgridBoundary(CellMask sgrd, SubGridBoundaryModes subGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge) {
                if (!object.ReferenceEquals(sgrd.GridData, this.GridData))
                    throw new ArgumentException("grid mismatch");
                m_SubGrid_InCells = sgrd;
                m_SubGridBoundaryTreatment = subGridBoundaryTreatment;
            }



            /// <summary>
            /// <see cref="CodomainMapping"/>
            /// </summary>
            UnsetteledCoordinateMapping m_CodomainMapping;

            /// <summary>
            /// coordinate mapping for the codomain variables;
            /// </summary>
            public UnsetteledCoordinateMapping CodomainMapping {
                get {
                    return m_CodomainMapping;
                }
            }

            /// <summary>
            /// <see cref="DomainMapping"/>
            /// </summary>
            UnsetteledCoordinateMapping m_DomainMapping;

            /// <summary>
            /// <see cref="Parameters"/>
            /// </summary>
            DGField[] m_Parameters;

            /// <summary>
            /// parameter mapping
            /// </summary>
            public IList<DGField> Parameters {
                get {
                    return m_Parameters;
                }
            }

            /// <summary>
            /// coordinate mapping for the domain variables;
            /// </summary>
            public UnsetteledCoordinateMapping DomainMapping {
                get {
                    return m_DomainMapping;
                }
            }

            /// <summary>
            /// Grid, on which this evaluator operates on.
            /// </summary>
            public IGridData GridData {
                get {
                    Debug.Assert(object.ReferenceEquals(DomainMapping.GridDat, CodomainMapping.GridDat));
                    return m_DomainMapping.GridDat;
                }
            }




            public double m_time = 0.0;

            /// <summary>
            /// Time passed e.g. to <see cref="CommonParams.time"/>, <see cref="CommonParamsBnd.time"/> and <see cref="CommonParamsVol.time"/>.
            /// </summary>
            public double time {
                get {
                    return m_time;
                }
                set {
                    m_time = value;
                }
            }

            /// <summary>
            /// Should return all DG fields which should be exchanged, see <see cref="m_TRX"/>.
            /// </summary>
            abstract protected DGField[] GetTrxFields();

            bool m_MPITtransceive = false;

            /// <summary>
            /// Turn MPI sending/receiving of parameters and domain fields on/off.
            /// </summary>
            public bool MPITtransceive {
                get {
                    //return m_TRX != null;
                    Debug.Assert((m_MPITtransceive == (m_TRX != null)) || (this.GetTrxFields().Length <= 0));

                    return m_MPITtransceive;
                }
                set {
                    m_MPITtransceive = value;

                    //if((m_TRX != null) && (value == true)) {
                    //    ArrayTools.ListEquals(m_TRX.)
                    //}


                    if((m_TRX == null) && (value == true) && (this.GetTrxFields().Length > 0)) {
                        // + + + + + + + + + +
                        // create transceiver
                        // + + + + + + + + + +

                        m_TRX = new Transceiver(GetTrxFields());
                    } else {
                        // + + + + + + + + + + + + +
                        // no communication required.
                        // + + + + + + + + + + + + +
                        m_TRX = null;
                    }

                }
            }



            /// <summary>
            /// Transceiver for the fields within <see cref="DomainMapping"/>
            /// </summary>
            protected Transceiver m_TRX;


        }



        /// <summary>
        /// evaluation of operators
        /// </summary>
        protected class EvaluatorNonLin : EvaluatorBase, IEvaluatorNonLin {

            /// <summary>
            /// Returns domain fields and parameters.
            /// </summary>
            protected override DGField[] GetTrxFields() {
                return ArrayTools.Cat(DomainFields.Fields, (base.Parameters != null) ? base.Parameters : new DGField[0]);
            }

            /// <summary>
            /// if the right-hand-side is present and contains nonlinear components, 
            /// this is the corresponding edge quadrature;
            /// otherwise, this member is null;
            /// </summary>
            BoSSS.Foundation.Quadrature.NonLin.NECQuadratureEdge m_NonlinearEdge;

            /// <summary>
            /// if the right-hand-side is present and contains nonlinear components, 
            /// this is the corresponding volume quadrature;
            /// otherwise, this member is null;
            /// </summary>
            BoSSS.Foundation.Quadrature.NonLin.NECQuadratureVolume m_NonlinearVolume;


            public EvaluatorNonLin(
                SpatialOperator owner,
                CoordinateMapping DomainVarMap,
                IList<DGField> ParameterMap,
                UnsetteledCoordinateMapping CodomainVarMap,
                EdgeQuadratureScheme edgeQrCtx,
                CellQuadratureScheme volQrCtx) //
             : base(owner, DomainVarMap, ParameterMap, CodomainVarMap) //
            {

                var grdDat = base.GridData;
                DomainFields = DomainVarMap;

                if(Owner.RequiresEdgeQuadrature) {


                    m_NonlinearEdge = new BoSSS.Foundation.Quadrature.NonLin.NECQuadratureEdge(grdDat,
                                                            Owner,
                                                            DomainVarMap,
                                                            ParameterMap,
                                                            CodomainVarMap,
                                                            edgeQrCtx.SaveCompile(grdDat, order));



                }

                if(Owner.RequiresVolumeQuadrature) {
                    m_NonlinearVolume = new BoSSS.Foundation.Quadrature.NonLin.NECQuadratureVolume(grdDat,
                                                                Owner,
                                                                DomainVarMap,
                                                                ParameterMap,
                                                                CodomainVarMap,
                                                                volQrCtx.SaveCompile(grdDat, order));

                }


                base.MPITtransceive = true;
            }

            /// <summary>
            /// DG filed which serve a input for the spatial operator.
            /// </summary>
            public CoordinateMapping DomainFields {
                get;
                private set;
            }
            
            /// <summary>
            /// evaluates the differential operator (<see cref="Owner"/>)
            /// for the domain variables/fields in <see cref="DomainMapping"/>, i.e.
            /// performs the operation
            /// <paramref name="output"/> = <paramref name="output"/>*<paramref name="beta"/> + Op(%)*<paramref name="alpha"/>
            /// </summary>
            /// <param name="output">
            /// (output): 
            /// The result of the operator evaluation,
            /// multiplied by <paramref name="alpha"/>,
            /// is <b>ACCUMULATED</b> here;
            /// It's up to the user to ensure that this array is initialized to 0.0,
            /// if necessary;
            /// Indices into this vector are computed according to <see cref="CodomainMapping"/>;
            /// </param>
            /// <param name="outputBndEdge">
            /// Some additional output vector for boundary fluxes, used only by the local time stepping
            /// </param>
            /// <param name="alpha">
            /// scaling of the operator;
            /// </param>
            /// <param name="beta">
            /// scaling applied to the accumulator; 
            /// </param>
            /// <remarks>
            /// This operation invokes MPI communication if <see cref="IEvaluator.MPITtransceive"/> is set to true: values of external cells are updated before
            /// fluxes are evaluated;
            /// </remarks>
            public void Evaluate<Tout>(double alpha, double beta, Tout output, double[] outputBndEdge = null)
                where Tout : IList<double> {
                using(var tr = new FuncTrace()) {
                    if(base.MPITtransceive)
                        MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                    if(output.Count != base.CodomainMapping.LocalLength)
                        throw new ArgumentOutOfRangeException("mismatch in length of output vector");

                    if(beta != 1.0)
                        GenericBlas.dscal(base.CodomainMapping.LocalLength, beta, output, 1);

                    base.UpdateCoefficients();

                    Debug.Assert((base.m_TRX != null) == (base.MPITtransceive == true));
                    if(base.m_TRX != null)
                        m_TRX.TransceiveStartImReturn();
#if DEBUG
                    output.CheckForNanOrInfV(true, true, true);
#endif

                    if(m_NonlinearVolume != null) {
                        using(new BlockTrace("Volume_Integration_NonLin", tr)) {
                            // volume integrals can be evaluated without knowing external cells
                            m_NonlinearVolume.m_Output = output;
                            m_NonlinearVolume.m_alpha = alpha;
                            m_NonlinearVolume.Time = time;
                            m_NonlinearVolume.Execute();
                            m_NonlinearVolume.m_Output = null;
                            m_NonlinearVolume.m_alpha = 1.0;
                        }

                    }



#if DEBUG
                    output.CheckForNanOrInfV(true, true, true);
#endif

                    if(base.m_TRX != null)
                        m_TRX.TransceiveFinish();

                    if(m_NonlinearEdge != null) {
                        using(new BlockTrace("Edge_Integration_NonLin", tr)) {
                            m_NonlinearEdge.m_Output = output;
                            m_NonlinearEdge.m_alpha = alpha;
                            m_NonlinearEdge.Time = time;
                            m_NonlinearEdge.SubGridBoundaryTreatment = base.SubGridBoundaryTreatment;
                            m_NonlinearEdge.SubGridCellsMarker = (base.m_SubGrid_InCells != null) ? base.m_SubGrid_InCells.GetBitMaskWithExternal() : null;

                            m_NonlinearEdge.m_outputBndEdge = outputBndEdge;

                            m_NonlinearEdge.Execute();

                            m_NonlinearEdge.m_Output = null;
                            m_NonlinearEdge.m_outputBndEdge = null;
                            m_NonlinearEdge.m_alpha = 1.0;
                            m_NonlinearEdge.SubGridCellsMarker = null;
                        }
                    }


#if DEBUG
                    output.CheckForNanOrInfV(true, true, true);
#endif

                }

            }


        }


        /// <summary>
        /// matrix assembly for linear or linearized operators
        /// </summary>
        protected class EvaluatorLinear : EvaluatorBase, IEvaluatorLinear {

            public EvaluatorLinear(
                SpatialOperator owner,
                UnsetteledCoordinateMapping DomainVarMap,
                IList<DGField> ParameterMap,
                UnsetteledCoordinateMapping CodomainVarMap,
                EdgeQuadratureScheme edgeQrCtx,
                CellQuadratureScheme volQrCtx) //
                 : base(owner, DomainVarMap, ParameterMap, CodomainVarMap) //
            {
                this.edgeRule = edgeQrCtx.SaveCompile(base.GridData, order);
                this.volRule = volQrCtx.SaveCompile(base.GridData, order);
                base.MPITtransceive = true;
            }

            /// <summary>
            /// Quadrature rule used for edge integration
            /// </summary>
            ICompositeQuadRule<QuadRule> edgeRule;

            /// <summary>
            /// quadrature rule used for volume integration
            /// </summary>
            ICompositeQuadRule<QuadRule> volRule;



            public void ComputeAffine<V>(V AffineOffset) where V : IList<double> {
                Internal_ComputeMatrixEx(default(BlockMsrMatrix), AffineOffset, true);
            }

            public void ComputeMatrix<M, V>(M Matrix, V AffineOffset)
                where M : IMutableMatrixEx
                where V : IList<double> // 
            {
                Internal_ComputeMatrixEx(Matrix, AffineOffset, false);
            }

            /// <summary>
            /// returns parameter fields
            /// </summary>
            protected override DGField[] GetTrxFields() {
                if(Parameters != null) {
                    return new DGField[0];
                } else {
                    List<DGField> FieldsForTransciever = new List<DGField>(Parameters.Count);
                    foreach(DGField ff in Parameters)
                        if(ff != null)
                            FieldsForTransciever.Add(ff);
                    return FieldsForTransciever.ToArray();
                }
            }

            /// <summary>
            /// matrix evaluation
            /// </summary>
            virtual protected void Internal_ComputeMatrixEx<M, V>(
                M Matrix, V AffineOffset, bool OnlyAffine)
                where M : IMutableMatrix
                where V : IList<double> //
            {
                MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
                using(var tr = new ilPSP.Tracing.FuncTrace()) {

                    // check arguments
                    // ===============
                    {

                        if(AffineOffset != null && AffineOffset.Count < base.CodomainMapping.LocalLength)
                            throw new ArgumentException("vector to short", "AffineOffset");
                        if(OnlyAffine == false) {
                            if(!Matrix.RowPartitioning.EqualsPartition(base.CodomainMapping))
                                throw new ArgumentException("wrong number of columns in matrix.", "Matrix");
                            if(!Matrix.ColPartition.EqualsPartition(base.DomainMapping))
                                throw new ArgumentException("wrong number of rows in matrix.", "Matrix");
                        }


                        if(base.m_SubGrid_InCells != null) {
                            throw new NotImplementedException("Subgrid feature is not implemented for matrix assembly");
                        }

                        base.UpdateCoefficients();

                    }

                    // do work
                    // =======

                    // transceiver is only needed for parameters
                    if(base.MPITtransceive) {
                        if(GetTrxFields().Length > 0) {
                            base.m_TRX.TransceiveStartImReturn();
                            base.m_TRX.TransceiveFinish();
                        }
                    }


                    // volume integration
                    // ------------------
                    if(volRule.Any()
                        && Owner.ContainesComponentType(typeof(IVolumeForm), typeof(IVolumeForm_UxV), typeof(IVolumeForm_UxGradV), typeof(IVolumeForm_GradUxV), typeof(IVolumeForm_GradUxGradV))) {
                        using(new BlockTrace("Volume_Integration_(new)", tr)) {
                            var mtxBuilder = new LECVolumeQuadrature2<M, V>(this.Owner);

                            mtxBuilder.Execute(volRule,
                                CodomainMapping, Parameters, DomainMapping,
                                OnlyAffine ? default(M) : Matrix, AffineOffset, time);

                        }

                    } else {
                        tr.Info("volume integration skipped: cell mask is empty");
                    }

                    // edge integration
                    // ----------------
                    if(!(edgeRule.IsNullOrEmpty())
                         && Owner.ContainesComponentType(typeof(IEdgeForm), typeof(IEdgeform_UxV), typeof(IEdgeform_UxGradV), typeof(IEdgeform_UxV), typeof(IEdgeSource_V))) {

                        using(new BlockTrace("Edge_Integration_(new)", tr)) {
                            var mxtbuilder2 = new LECEdgeQuadrature2<M, V>(this.Owner);

                            mxtbuilder2.Execute(edgeRule, CodomainMapping, Parameters, DomainMapping, OnlyAffine ? default(M) : Matrix, AffineOffset, time);
                            tr.Info("done.");
                            mxtbuilder2 = null;
                        }
                    }
                }
            }
        }
    }
}
