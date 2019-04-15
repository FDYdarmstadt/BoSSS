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
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// An operator which is specialized in XDG Fields, i.e.
    /// it can have components which couple the phases.
    /// Mk2: enables the definition of different equationcomponents for each phase 
    /// </summary>
    public class XSpatialOperatorMk2 {

        /// <summary>
        /// Non-coupling surface terms; originally intended to implement the flux-form of the surface tension.
        /// </summary>
        public SpatialOperator SurfaceElementOperator {
            get;
            private set;
        }


        //void ConstructorCommon() {
        //    SurfaceElementOperator = new FixedOrder_SpatialOperator(base.DomainVar, base.ParameterVar, base.CodomainVar);
        //}


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
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(int,int,int,Func{int[],int[],int[],int},string[])"/>
        /// </summary>
        public XSpatialOperatorMk2(int NoOfDomFields, int NoOfParameters, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, params string[] __varnames)
            : this(GetSubarray(__varnames, 0, NoOfDomFields), GetSubarray(__varnames, NoOfDomFields, NoOfParameters), GetSubarray(__varnames, NoOfDomFields + NoOfParameters, NoOfCodomFields), QuadOrderFunc) {
            if(NoOfCodomFields + NoOfDomFields + NoOfParameters != __varnames.Length)
                throw new ArgumentException("mismatch between number of provided variable names and given number of domain, parameter and codomain fields.");

        }

        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(int,int,Func{int[],int[],int[],int},string[])"/>
        /// </summary>
        public XSpatialOperatorMk2(int NoOfDomFields, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, params string[] __varnames)
           : this(GetSubarray(__varnames, 0, NoOfDomFields), GetSubarray(__varnames, NoOfDomFields, NoOfCodomFields), QuadOrderFunc) {
            if(NoOfCodomFields + NoOfDomFields != __varnames.Length)
                throw new ArgumentException("mismatch between number of provided variable names and given number of domain and codomain fields.");
        }

        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(IList{string},IList{string},Func{int[],int[],int[],int})"/>
        /// </summary>
        public XSpatialOperatorMk2(IList<string> __DomainVar, IList<string> __CoDomainVar, Func<int[], int[], int[], int> QuadOrderFunc)
            : this(__DomainVar, null, __CoDomainVar, QuadOrderFunc) {
        }

        /// <summary>
        /// ctor, see <see cref="SpatialOperator.SpatialOperator(IList{string},IList{string},IList{string},Func{int[],int[],int[],int})"/>
        /// </summary>
        public XSpatialOperatorMk2(IList<string> __DomainVar, IList<string> __ParameterVar, IList<string> __CoDomainVar, Func<int[], int[], int[], int> QuadOrderFunc) {
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
            m_EquationComponentsHelper = new _XEquationComponents(this);
            this.QuadOrderFunction = QuadOrderFunc;


            SurfaceElementOperator = new FixedOrder_SpatialOperator(DomainVar, ParameterVar, CodomainVar);

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


        /// <summary>
        /// <see cref="EquationComponents"/>
        /// </summary>
        SortedList<string, List<IEquationComponent>> m_EquationComonents;

        _XEquationComponents m_EquationComponentsHelper;

        /// <summary>
        /// for each variable in <see cref="CodomainVar"/>, a
        /// collection of equation components that define the operator.
        /// 
        /// </summary>
        public _XEquationComponents EquationComponents {
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

            SurfaceElementOperator.Commit();
        }


        /// <summary>
        /// implementation of <see cref="EquationComponents" />;
        /// </summary>
        public class _XEquationComponents : IEnumerable<KeyValuePair<string, IEnumerable<IEquationComponent>>> {

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




        class FixedOrder_SpatialOperator : SpatialOperator {
            public FixedOrder_SpatialOperator(IList<string> __DomainVar, IList<string> __ParameterVar, IList<string> __CoDomainVar)
                : base(__DomainVar, __ParameterVar, __CoDomainVar, null) //
            {
                base.QuadOrderFunction = this.QOF;
            }

            int QOF(int[] A, int[] B, int[] C) {
                return m_Order;
            }

            internal int m_Order;
        }

    }

}
