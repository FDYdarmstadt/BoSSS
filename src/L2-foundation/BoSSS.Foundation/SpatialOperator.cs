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
    /// Delegate to trigger the update of parameter fields (e.g. when computing finite difference Jacobian, see e.g. <see cref="SpatialOperator.GetFDJacobianBuilder(IList{DGField}, IList{DGField}, UnsetteledCoordinateMapping, EdgeQuadratureScheme, CellQuadratureScheme)"/>).
    /// </summary>
    /// <param name="DomainVar">
    /// Input fields, current state of domain variables
    /// </param>
    /// <param name="ParameterVar">
    /// Output fields, updated states of parameter fields
    /// </param>
    public delegate void DelParameterUpdate(IEnumerable<DGField> DomainVar, IEnumerable<DGField> ParameterVar);

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
        internal protected void Verify() {
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
        /// are concatenated with domain variable names (see <see cref="IEquationComponent.ArgumentOrdering"/>).
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
        /// constructs a new evaluator object for explicit evaluation this spatial operator
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

        /// <summary>
        /// Creator of a <see cref="EvaluatorLinear"/> object.
        /// </summary>
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

                    if (!object.ReferenceEquals(DomainVarMap.GridDat, CodomainVarMap.GridDat))
                        throw new ArgumentException("Domain and Codomain map must be assigned to the same grid");

                    foreach(var f in Parameters) {
                        if(f != null) {
                            if (!object.ReferenceEquals(DomainVarMap.GridDat, f.GridDat))
                                throw new ArgumentException("Parameter fields, domain and codomain basis must be assigned to the same grid");
                        }
                    }


                    m_Owner = owner;
                    m_CodomainMapping = CodomainVarMap;
                    m_DomainMapping = DomainVarMap;
                    m_Parameters = (ParameterMap != null) ? ParameterMap.ToArray() : new DGField[0];


                    //IEnumerable<Basis> allBasis = DomainVarMap.BasisS;
                    //if (ParameterMap != null) {
                    //    allBasis = allBasis.Union(ParameterMap.Select(f => f.Basis));
                    //}
                    //allBasis = allBasis.Union(CodomainVarMap.BasisS);
                    //IGridData grdDat = allBasis.First().GridDat;
                    //foreach (var b in allBasis) {
                    //    if (!object.ReferenceEquals(grdDat, b.GridDat)) {
                    //        throw new ArgumentException("all fields (domain, parameter, codomain) must be defined on the same grid.");
                    //    }
                    //}

                    if (!m_Owner.IsCommited)
                        throw new ApplicationException("operator assembly must be finalized before by calling 'Commit' before this method can be called.");

                    order = owner.GetOrderFromQuadOrderFunction(m_DomainMapping, ParameterMap, CodomainVarMap);

                    m_OperatorCoefficients = new CoefficientSet() {
                        CellLengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Cells.CellLengthScale,
                        EdgeLengthScales = ((BoSSS.Foundation.Grid.Classic.GridData)(this.GridData)).Edges.h_min_Edge,
                        UserDefinedValues = new Dictionary<string, object>(),
                        GrdDat = this.GridData
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
                    var RealTrxFields = this.GetTrxFields().Where(f => f != null).ToArray();

                    //Debug.Assert((m_MPITtransceive == (m_TRX != null)) || (RealTrxFields.Length <= 0));

                    return m_MPITtransceive;
                }
                set {
                    m_MPITtransceive = value;

                    //if((m_TRX != null) && (value == true)) {
                    //    ArrayTools.ListEquals(m_TRX.)
                    //}
                    var RealTrxFields = this.GetTrxFields().Where(f => f != null).ToArray();
                    
                    if((value == true) && (RealTrxFields.Length > 0)) {
                        // + + + + + + + + + +
                        // create transceiver
                        // + + + + + + + + + +

                        if(m_TRX == null)
                            m_TRX = new Transceiver(RealTrxFields);
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

            /// <summary>
            /// Not for direct user interaction
            /// </summary>
            internal protected EvaluatorNonLin(
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

                    //if(!rem) {
                    //    rem = true;
                    //    Console.WriteLine("Reminder: edge terms deactivated.");
                    //}

                    
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

            /// <summary>
            /// Not for direct user interaction
            /// </summary>
            internal protected EvaluatorLinear(
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


            /// <summary>
            /// For the operator linearization 
            /// \f[
            ///    \mathcal{M} U + \mathcal{B}
            /// \f]
            /// this computes only the vector \f$ \mathcal{B} \f$
            /// </summary>
            /// <param name="AffineOffset"></param>
            public void ComputeAffine<V>(V AffineOffset) where V : IList<double> {
                Internal_ComputeMatrixEx(default(BlockMsrMatrix), AffineOffset, true);
            }

            /// <summary>
            /// computes a linearization of the operator in the form 
            /// \f[
            ///    \mathcal{M} U + \mathcal{B}.
            /// \f]
            /// </summary>
            /// <param name="Matrix">
            /// Output, the operator matrix is accumulated here
            /// </param>
            /// <param name="AffineOffset">
            /// Output, the affine part of the operator linearization is accumulated here
            /// </param>
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
                            mxtbuilder2 = null;
                            //Console.WriteLine("edge lin deact");
                        }
                    }
                }
            }
        }



        /// <summary>
        /// constructs a <see cref="FDJacobianBuilder"/> object to linearize nonlinear operators
        /// </summary>
        public virtual IEvaluatorLinear GetFDJacobianBuilder(
            IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
            DelParameterUpdate __delParameterUpdate,
            EdgeQuadratureScheme edgeQrCtx = null,
            CellQuadratureScheme volQrCtx = null) //
        {

            using(new FuncTrace()) {

                /// This is already done in the constructor of Evaluator
#if DEBUG
                if(!m_IsCommited)
                    throw new ApplicationException("operator assembly must be finalized before by calling 'Commit' before this method can be called.");
#endif


                var e = new FDJacobianBuilder(new EvaluatorNonLin(
                    this, 
                    new CoordinateMapping(DomainFields), ParameterMap, CodomainVarMap, 
                    edgeQrCtx, volQrCtx),
                    __delParameterUpdate);
                    //new CoordinateMapping(DomainFields), ParameterMap, CodomainVarMap, edgeQrCtx, volQrCtx);

                return e;
            }
        }


        /// <summary>
        /// Computes the (approximate) Jacobian matrix of the spatial operator by finite differences.
        /// </summary>
        public class FDJacobianBuilder : IEvaluatorLinear {

            /// <summary>
            /// Not for direct user interaction
            /// </summary>
            internal protected FDJacobianBuilder(IEvaluatorNonLin __Eval, DelParameterUpdate __delParameterUpdate) {

                eps = 1.0;
                while( 1.0 + eps > 1.0) {
                    eps = eps /2;
                }
                eps = Math.Sqrt(eps);

                Eval = __Eval;
                DelParamUpdate = __delParameterUpdate;

                BuildGridColoring();
            }

            /// <summary>
            /// Approximately the square root of the machine epsilon.
            /// </summary>
            double eps;

            /// <summary>
            /// Epsilon used for the finite difference
            /// </summary>
            public double Eps {
                get {
                    return eps;
                }
                set {
                    if(value <= 0.0)
                        throw new ArgumentOutOfRangeException();
                    eps = value;
                }
            }

            /// <summary>
            /// Grid
            /// </summary>
            public IGridData GridData {
                get {
                    return Eval.GridData;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public UnsetteledCoordinateMapping CodomainMapping {
                get {
                    return Eval.CodomainMapping;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public UnsetteledCoordinateMapping DomainMapping {
                get {
                    return Eval.DomainMapping;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public IList<DGField> Parameters {
                get {
                    return Eval.Parameters;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public SubGridBoundaryModes SubGridBoundaryTreatment {
                get {
                    return Eval.SubGridBoundaryTreatment;
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public double time {
                get { return Eval.time; }
                set { Eval.time = value; }
            }

            /// <summary>
            /// 
            /// </summary>
            public SpatialOperator Owner {
                get {
                    return Eval.Owner;
                }
            }

            /// <summary>
            /// can be dangerous to turn off
            /// </summary>
            public bool MPITtransceive {
                get => throw new NotImplementedException();
                set => throw new NotImplementedException();
            }

            /// <summary>
            /// 
            /// </summary>
            virtual public CoefficientSet OperatorCoefficients {
                get { return Eval.OperatorCoefficients; }
                set { OperatorCoefficients = value; }
            }

            IEvaluatorNonLin Eval;


            DelParameterUpdate DelParamUpdate;

            /// <summary>
            /// - 1st index: enumeration of color lists
            /// - 2nd index: enumeration of cells in color list
            /// - content: local cell index in which an Epsilon-distortion is applied
            /// </summary>
            int[][] ColorLists;

            /// <summary>
            /// Cells which are distorted on an other processor, but influence the result on this processor
            /// - 1st index: correlates to 1st index of <see cref="ColorLists"/>
            /// - 2nd index: enumeration of cells
            /// - content: some external cell index; this determines the column index of the finite difference result into the Jacobian matrix
            /// </summary>
            int[][] ExternalColorLists;
            
            /// <summary>
            /// - 1st index: correlates to 1st index of <see cref="ColorLists"/>
            /// - 2nd index: correlates with 2nd index of <see cref="ExternalColorLists"/>
            /// - 3rd index: enumeration
            /// - content: a local cell index of some cell which is affected by an epsilon-distortion on some other processor; this determines the row index of the finite difference result into the Jacobian matrix
            /// </summary>
            int[][][] ExternalColorListsNeighbors;


            void BuildGridColoring() {
                var gDat = Eval.GridData;


                int[][] Neighs = gDat.iLogicalCells.CellNeighbours;
                int J = gDat.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = gDat.iLogicalCells.Count;
                long[] GlidxExt = gDat.iParallel.GlobalIndicesExternalCells;
                var Gl2LocExt = gDat.iParallel.Global2LocalIdx;
                var CellPart = gDat.CellPartitioning;

                Random rnd = new Random(gDat.MpiRank + 1);

#if DEBUG
                for (int j = 0; j < J; j++) {
                    int[] CN = Neighs[j];
                    foreach(int jN in CN) {
                        Debug.Assert(jN >= 0);
                        Debug.Assert(jN != j);
                    }
                }
#endif

                int[] LocalMarker = new int[JE]; //    marker for blocked in the current pass 
                int[] ExchangedMarker = new int[JE]; //  accumulation buffer for MPI exchange
                BitArray Colored = new BitArray(JE); // all cells which are already colored (in previous passes)
                BitArray ColoredPass = new BitArray(JE); // all cells which are colored in current pass
                int[] LocalColorCause = new int[JE];

                List<int> CellList = new List<int>();
                List<int[]> ColorListsTmp = new List<int[]>();
                List<int[]> ExternalColorListsTmp = new List<int[]>();
                List<int[][]> ExternalColorListsNeighborsTmp = new List<int[][]>();

                int myMarkerToken = gDat.MpiRank + 1;
                int bRun = 0xFFFFFF;
                int locColoredCells = 0;
                int DeadlockWatch = 0;
                while(bRun != 0) {

                    // find next color list
                    // ====================

                    Array.Clear(LocalMarker, 0, LocalMarker.Length);// LocalMarker.SetAll(int.MaxValue);
                    ColoredPass.SetAll(false);

                    for(int j = 0; j < J; j++) {
                        if (Colored[j] == true)
                            continue;
                        if (LocalMarker[j] != 0)
                            continue;

                        int[] Neighs_j = Neighs[j];
                        Debug.Assert(Neighs_j.Contains(j) == false, "Cell seems to be its own neighbor.");
                        bool cont = false;
                        foreach(int jn in Neighs_j) {
                            if(LocalMarker[jn] != 0) {
                                cont = true;
                                break;
                            }
                        }
                        if(cont)
                            continue;

                        // if we reached this point, we finally found a cell which we are allowed to add to the current color set.
                        ColoredPass[j] = true;                        
                        LocalMarker[j] = myMarkerToken;
                        LocalColorCause[j] = j;
                        foreach(int jn in Neighs_j) {
                            LocalMarker[jn] = -myMarkerToken; // mark neighbor cells with negative token!
                            LocalColorCause[jn] = j;
                        }
                    }

                    // fix parallel conflicts
                    // ======================

                    var LocalMarker_Bkup = LocalMarker.CloneAs();
                    var Removed = new List<int>();
                    int[] ExchangedMarker_Bkup = null;

                    if(gDat.MpiSize > 1) {
                        ExchangedMarker_Bkup = ExchangedMarker.CloneAs();

                        int GlobalConflicts = 999;
                        do {
                            Array.Copy(LocalMarker, ExchangedMarker, JE);
                            ExchangedMarker.MPIExchange(gDat);

                            int LocalConflicts = 0;

                            for (int je = J; je < JE; je++) {
                                if (LocalMarker[je] != 0 && ExchangedMarker[je] != 0) {
                                    Debug.Assert(LocalMarker[je] != ExchangedMarker[je]);
                                    LocalConflicts++;

                                    double rndVal = rnd.NextDouble();

                                    //// some parallel conflict detected: one of the two ranks has to yield

                                    //if (ExchangedMarker[je] > 0 && Math.Abs(ExchangedMarker[je]) > myMarkerToken) {
                                    //    // the other rank should yield
                                    //} else {
                                    //    // this rank has to yield
                                    if (rndVal >= 0.5) {
                                        int jToRemove = LocalColorCause[je];
                                        Debug.Assert(jToRemove < J);
                                        Debug.Assert(ColoredPass[jToRemove] == true);

                                        Removed.Add(jToRemove);

                                        ColoredPass[jToRemove] = false;
                                        LocalMarker[jToRemove] = 0;
                                        int[] Neighs_jToRemove = Neighs[jToRemove];
                                        foreach (int jn in Neighs_jToRemove) {
                                            LocalMarker[jn] = 0;
                                        }

                                        //}
                                    }
                                }
                            }

                            GlobalConflicts = LocalConflicts.MPISum();

                        } while (GlobalConflicts > 0);

#if DEBUG
                        Array.Copy(LocalMarker, ExchangedMarker, JE);
                        ExchangedMarker.MPIExchange(gDat);
                        for(int je = J; je < JE; je++) {
                            Debug.Assert(ExchangedMarker[je] == 0 || LocalMarker[je] == 0, "Error in conflict resolution.");
                        }
#endif
                    }



                    // remember recently found color list
                    // ==================================
                    CellList.Clear();
                    for(int j = 0; j < J; j++) {
                        if(ColoredPass[j]) {
                            locColoredCells++;
                            CellList.Add(j);
                            Debug.Assert(Colored[j] == false);
                            Colored[j] = true;
                        }
                    }
                    int LocColoredPass = CellList.Count;

                    int GlobColoredPass = LocColoredPass.MPISum();
                    //Console.WriteLine("Colored in pass {0}: {1}", ColorListsTmp.Count, GlobColoredPass);
                    if (GlobColoredPass <= 0) {
                        DeadlockWatch++;
                        if(DeadlockWatch >= 1000)
                            throw new ApplicationException("Deadlock in parallel coloring.");
                        continue;
                        //Debugger.Launch();
                    }


                    ColorListsTmp.Add(CellList.ToArray());

                    // communicate external lists
                    // ==========================

                    if (gDat.MpiSize > 1) {

                        //Debugger.Launch();

                        var ExchData = new Dictionary<int, List<Tuple<int, int>>>();

                        foreach(int j in CellList) {
                            int[] Neighs_j = Neighs[j];
                            foreach (int jN in Neighs_j) {
                                if(jN >= J) {
                                    
                                    int Gl_jN = (int) GlidxExt[jN - J];
                                    int iProc = CellPart.FindProcess(Gl_jN);
                                    int Gl_j = j + CellPart.i0;

                                    if(!ExchData.TryGetValue(iProc, out var ExchData_iProc)) {
                                        ExchData_iProc = new List<Tuple<int, int>>();
                                        ExchData.Add(iProc, ExchData_iProc);
                                    }

                                    ExchData_iProc.Add(new Tuple<int, int>(Gl_j, Gl_jN));
                                }
                            }
                        }

                        var RcvData = SerialisationMessenger.ExchangeData(ExchData);

                        var ExtColor = new Dictionary<int, List<int>>();

                        foreach (var kv in RcvData) {
                            int iProc = kv.Key;
                            var list = kv.Value;

                            foreach(var t in list) {
                                int Gl_j = t.Item1;
                                int Gl_jN = t.Item2;
                                Debug.Assert(CellPart.FindProcess(Gl_j) == iProc);
                                Debug.Assert(CellPart.IsInLocalRange(Gl_jN));

                                int Loc_jN = Gl_jN - CellPart.i0;
                                Debug.Assert(Loc_jN >= 0 && Loc_jN < J);
                                int Loc_j = Gl2LocExt[Gl_j];
                                Debug.Assert(Loc_j >= J && Loc_j < JE);

                                List<int> Neighs_Loc_jN;
                                if(!ExtColor.TryGetValue(Loc_j, out Neighs_Loc_jN)) {
                                    Neighs_Loc_jN = new List<int>();
                                    ExtColor.Add(Loc_j, Neighs_Loc_jN);
                                }

                                Neighs_Loc_jN.Add(Loc_jN);
                            }
                        }


                        int[] T2 = new int[ExtColor.Count];
                        int[][] T3 = new int[ExtColor.Count][];
                        int cnt = 0;
                        foreach(var kv in ExtColor) {
                            T2[cnt] = kv.Key;
                            T3[cnt] = kv.Value.ToArray();
                            cnt++;
                        }

                        ExternalColorListsTmp.Add(T2);
                        ExternalColorListsNeighborsTmp.Add(T3);

                    } else {
                        ExternalColorListsTmp.Add(new int[0]);
                        ExternalColorListsNeighborsTmp.Add(new int[0][]);
                    }


                    // check for loop termination
                    // ==========================
                    Debug.Assert(locColoredCells <= J);
                    int bRunLoc = 0xFFFFFF;
                    if(locColoredCells >= J) 
                        bRunLoc = 0;
                    bRun = bRunLoc.MPIMax();
                }

                // store
                // =====
                this.ColorLists = ColorListsTmp.ToArray();
                this.ExternalColorLists = ExternalColorListsTmp.ToArray();
                this.ExternalColorListsNeighbors = ExternalColorListsNeighborsTmp.ToArray();

                // checks
                // ======
#if DEBUG
                {
                    int[] TouchLoc = new int[JE];

                    foreach(int[] _CellList in this.ColorLists) {

                        int[] NeighTouchLoc = new int[JE];

                        foreach(int j in _CellList) {
                            Debug.Assert(j < JE);
                            Debug.Assert(TouchLoc[j] == 0);
                            TouchLoc[j]++;

                            Debug.Assert(NeighTouchLoc[j] == 0);
                            NeighTouchLoc[j]++;

                            foreach(int jn in Neighs[j]) {
                                Debug.Assert(NeighTouchLoc[jn] == 0);
                                NeighTouchLoc[jn]++;
                            }

                        }

                        int[] NeighTouchGlob = NeighTouchLoc.CloneAs();
                        NeighTouchGlob.MPIExchange(gDat);
                        for(int j = J; j < JE; j++) {
                            NeighTouchGlob[j] += NeighTouchLoc[j];
                        }

                        for(int j = 0; j < JE; j++) {
                            Debug.Assert(NeighTouchGlob[j] <= 1);
                        }
                    }

                    int[] TouchGlob = TouchLoc.CloneAs();
                    TouchGlob.MPIExchange(gDat);
                    for(int j = J; j < JE; j++) {
                        TouchGlob[j] += TouchLoc[j];
                    }
                    
                    for(int j = 0; j < JE; j++) {
                        Debug.Assert(TouchGlob[j] == 1);
                    }

                }
#endif

            }


            /// <summary>
            /// computes a approximate linearization of the operator in the form 
            /// \f[
            ///    \mathcal{M} U + \mathcal{B}.
            /// \f]
            /// </summary>
            /// <param name="Matrix">
            /// Output, the approximate Jacobian matrix of the operator is accumulated here
            /// </param>
            /// <param name="AffineOffset">
            /// Output, the operator value in the linearization point
            /// </param>
            public void ComputeMatrix<M, V>(M Matrix, V AffineOffset)
                where M : IMutableMatrixEx
                where V : IList<double> // 
            {
                
                // init locals
                // ===========
                var codMap = Eval.CodomainMapping;
                var domMap = Eval.DomainMapping;
                DGField[] domFields = Eval.DomainFields.Fields.ToArray();
                var U0 = new CoordinateVector(Eval.DomainFields);

                int j0 = Eval.GridData.CellPartitioning.i0;
                int J = Eval.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = Eval.GridData.iLogicalCells.Count;
                int NoOfDomFields = domMap.BasisS.Count;
                int NoOfCodFields = codMap.BasisS.Count;

                int Lout = Eval.CodomainMapping.LocalLength;
                int Lin = domMap.LocalLength;

                int[][] Neighs = Eval.GridData.iLogicalCells.CellNeighbours;

                var lastCodB = codMap.BasisS.Last();
                var lastDomB = domMap.BasisS.Last();

                int NoOfEvals = 0;

                // check Args
                // ==========

                if(!Matrix.RowPartitioning.EqualsPartition(codMap))
                    throw new ArgumentException("Mismatch in matrix row partition.");
                if(!Matrix.ColPartition.EqualsPartition(domMap))
                    throw new ArgumentException("Mismatch in matrix column partition.");
                if(AffineOffset.Count != codMap.LocalLength)
                    throw new ArgumentException("Mismatch in length of affine offset.");

                // evaluate at linearization point
                // ===============================

                double[] F0 = new double[Lout];
                DelParamUpdate(domFields, Eval.Parameters.ToArray());
                Eval.Evaluate(1.0, 0.0, F0);
                AffineOffset.AccV(1.0, F0);
                NoOfEvals++;

                // compute epsilon's
                // =================

                double[] Epsilons = new double[domMap.Ntotal];
                double relEps = this.Eps;
                //double absEps = 1.0e-15; 
                double absEps = 1.0;
                for(int i = 0; i < Lin; i++) {
                    double EpsBase = Math.Abs(U0[i]);
                    if(EpsBase < absEps)
                        EpsBase = absEps;

                    Epsilons[i] = EpsBase * relEps;
                }
                Epsilons.MPIExchange(Eval.GridData);


                // compute directional derivatives
                // ===============================

                double[] U0backup = new double[Lin];
                double[] EvalBuf = new double[Lout];
                if(!domMap.AllBlockSizesEqual)
                    throw new NotSupportedException();
                MultidimensionalArray Buffer = MultidimensionalArray.Create(Lout, domMap.GetBlockLen(domMap.FirstBlock));

                for (int iCellPass = 0; iCellPass < ColorLists.Length; iCellPass++) { // loop over all cell lists...
                    int[] CellList = this.ColorLists[iCellPass];
                    int[] ExtCellList = this.ExternalColorLists[iCellPass];

                    int[] CoordCounter = new int[JE];
                    int[] FieldCounter = new int[JE];

                    int maxNj = 0;
                    foreach (int j in CellList) {
                        int Nj = domMap.GetTotalNoOfCoordinatesPerCell(j);
                        maxNj = Math.Max(Nj, maxNj);
                    }
                    maxNj = maxNj.MPIMax();

                    Buffer.Clear();

                    for (int n = 0; n < maxNj; n++) { // loop over DG coordinates in cell

                        // backup DG coordinates
                        // ---------------------
                        U0backup.SetV(U0);

                        // apply distortions
                        // -----------------
                        int AnyLoc = 0;
                        foreach (int j in CellList) {
                            int iFld = FieldCounter[j];
                            int nFld = CoordCounter[j];
                            if (iFld > NoOfDomFields)
                                continue; // finished with cell 'j'

                            AnyLoc = -1;

                            Debug.Assert(j >= 0 && j < J);
                            int iLoc = domMap.LocalUniqueCoordinateIndex(iFld, j, nFld);
                            Debug.Assert(iLoc >= 0 && iLoc < domMap.LocalLength);

                            double oldVal = domFields[iFld].Coordinates[j, nFld];
                            domFields[iFld].Coordinates[j, nFld] = oldVal + Epsilons[iLoc];
                            Debug.Assert(domFields[iFld].Coordinates[j, nFld] != oldVal);
                        }
                        int AnyGlob = AnyLoc.MPIMin();
                        if (AnyGlob >= 0)
                            break; // finished with entire cell list on all processors

                        // evaluate operator
                        // -------------------
                        EvalBuf.ClearEntries();
                        DelParamUpdate(domFields, Eval.Parameters.ToArray());
                        Eval.Evaluate(1.0, 0.0, EvalBuf);
                        NoOfEvals++;

                        // ------------------------------

                        for (int IntExt = 0; IntExt < 2; IntExt++) {
                            int[] __CellList;
                            switch (IntExt) {
                                case 0: __CellList = CellList; break;
                                case 1: __CellList = ExtCellList; break;
                                default: throw new ApplicationException();
                            }


                            // save results
                            // -------------------------------
                            int cnt = 0;
                            foreach (int _j in __CellList) {
                                int[] Neighs_j; // = Neighs[_j];
                                switch (IntExt) {
                                    case 0: Neighs_j = Neighs[_j]; break;
                                    case 1: Neighs_j = this.ExternalColorListsNeighbors[iCellPass][cnt]; break;
                                    default: throw new ApplicationException();
                                }
                                cnt++;

                                int jCol = _j;

                                int iFldCol = FieldCounter[jCol];
                                int nFldCol = CoordCounter[jCol];
                                if (iFldCol > NoOfDomFields)
                                    continue; // finished with cell

                                int iCol = domMap.LocalUniqueCoordinateIndex(iFldCol, jCol, nFldCol);
                                int i0Col = domMap.LocalUniqueCoordinateIndex(0, jCol, 0);
                                int iRelCol = iCol - i0Col;

                                for (int k = 0; k <= Neighs_j.Length; k++) { // loop over neighbors which are influenced by the distortion
                                    int jRow;
                                    if (k == 0) {
                                        jRow = _j;
                                    } else {
                                        jRow = Neighs_j[k - 1];
                                    }

                                    if (jRow >= J)
                                        continue; // external cell; should be treated on other proc.

                                    int i0Row = codMap.LocalUniqueCoordinateIndex(0, jRow, 0);
                                    int NoOfRows = codMap.GetBlockLen(jRow);

                                    for (int iRelRow = 0; iRelRow < NoOfRows; iRelRow++) {
                                        int iRow = i0Row + iRelRow;

                                        double u1 = EvalBuf[iRow];
                                        double u0 = F0[iRow];
                                        double h = Epsilons[iCol];

                                        double diff = (u1 - u0) / h;
                                        Buffer[iRow, iRelCol] = diff;
                                    }
                                }
                            }

                            // increase counters
                            // ------------------
                            foreach (int j in __CellList) {
                                int iFld = FieldCounter[j];
                                if (iFld > NoOfDomFields)
                                    continue; // finished with cell 'j'

                                int Nj = domMap.BasisS[iFld].GetLength(j);
                                CoordCounter[j]++;
                                if (CoordCounter[j] >= Nj) {
                                    CoordCounter[j] = 0;
                                    FieldCounter[j]++;
                                }

                            }
                        }

                        // restore original DG coordinates
                        // -------------------------------
                        U0.SetV(U0backup);
                    }

                    // save to matrix
                    // --------------

                    for (int IntExt = 0; IntExt < 2; IntExt++) {
                        int[] __CellList;
                        switch (IntExt) {
                            case 0: __CellList = CellList; break;
                            case 1: __CellList = ExtCellList; break;
                            default: throw new ApplicationException();
                        }

                        int cnt = 0;
                        foreach (int _j in __CellList) {
                            int[] Neighs_j; // = Neighs[_j];
                            switch (IntExt) {
                                case 0: Neighs_j = Neighs[_j]; break;
                                case 1: Neighs_j = this.ExternalColorListsNeighbors[iCellPass][cnt]; break;
                                default: throw new ApplicationException();
                            }
                            cnt++;

                            int jCol = _j;
                            int i0Col = domMap.LocalUniqueCoordinateIndex(0, jCol, 0);
                            int iECol = domMap.LocalUniqueCoordinateIndex(NoOfDomFields - 1, jCol, lastDomB.GetLength(jCol) - 1);

                            for (int k = 0; k <= Neighs_j.Length; k++) { // loop over neighbors which are influenced by the distortion
                                int jRow;
                                if (k == 0) {
                                    jRow = _j;
                                } else {
                                    jRow = Neighs_j[k - 1];
                                }

                                if (jRow >= J)
                                    continue; // external cell; should be treated on other proc.


                                int i0Row = domMap.LocalUniqueCoordinateIndex(0, jRow, 0);
                                int iERow = domMap.LocalUniqueCoordinateIndex(NoOfCodFields - 1, jRow, lastCodB.GetLength(jRow) - 1);

                                var Block = Buffer.ExtractSubArrayShallow(new int[] { i0Row, 0 }, new int[] { iERow, iECol - i0Col });

                                Matrix.AccBlock(i0Row + codMap.i0,
                                    //i0Col + domMap.i0, 
                                    domMap.GlobalUniqueCoordinateIndex(0, jCol, 0),
                                    1.0, Block);
                            }
                        }

                    }
                }
                // restore original state before return
                // ====================================
                U0.SetV(U0backup);
                DelParamUpdate(domFields, Eval.Parameters.ToArray());
                //Console.WriteLine("Total number of evaluations: " + NoOfEvals);
            }

            /// <summary>
            /// Evaluation at the linearization point
            /// </summary>
            public void ComputeAffine<V>(V AffineOffset) where V : IList<double> {
                int Lout = Eval.CodomainMapping.LocalLength;
               
                double[] F0 = new double[Lout];
                DelParamUpdate(Eval.DomainFields.Fields.ToArray(), Eval.Parameters.ToArray());
                Eval.Evaluate(1.0, 0.0, F0);
                AffineOffset.AccV(1.0, F0);
               
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="sgrd"></param>
            /// <param name="subGridBoundaryTreatment"></param>
            public void ActivateSubgridBoundary(CellMask sgrd, SubGridBoundaryModes subGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge) {
                Eval.ActivateSubgridBoundary(sgrd, subGridBoundaryTreatment);
            }
        }



    }
}
