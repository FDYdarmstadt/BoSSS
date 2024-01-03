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
    /// This class represents a spatial operator which maps
    /// from (DG-) variables in the domain, identified and ordered by <see cref="DomainVar"/>
    /// to variables in the co-domain, identified and ordered by <see cref="CodomainVar"/>.
    /// </summary>
    public class DifferentialOperator : IDifferentialOperator {


        bool m_IsLinear;

        /// <summary>
        /// true, if the PDE defined by operator can entirely be solved by a linear solver
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

                return Quadrature.PeriodicBoundaryUtils.GetVectorFieldIndices(this.DomainVar, 3);
            }
        }


        /// <summary>
        /// <see cref="IDifferentialOperator.SolverSafeguard"/>
        /// </summary>
        public SolverSafeguard SolverSafeguard {
            get;
            set;
        }

        Func<int[], int[], int[], int> m_QuadOrderFunction;

        /// <summary>
        /// Function Mapping from Domain Variable Degrees, Parameter Degrees and CoDomain Variable Degrees to the Quadrature Order
        /// </summary>
        public Func<int[], int[], int[], int> QuadOrderFunction {
            get {
                return m_QuadOrderFunction;
            }
            set {
                // deactivated due to legacy code issues:
                //if(IsCommited)
                //    throw new NotSupportedException("not allowed to change after Commit");
                m_QuadOrderFunction = value;
            }
        }

        Func<IGridData, EdgeQuadratureScheme> m_EdgeQuadraturSchemeProvider;

        /// <summary>
        /// User-customizable factory, to specify the edge quadrature, see also <see cref="QuadOrderFunction"/>
        /// </summary>
        public Func<IGridData,EdgeQuadratureScheme> EdgeQuadraturSchemeProvider {
            get {
                if(m_EdgeQuadraturSchemeProvider == null)
                    m_EdgeQuadraturSchemeProvider = (IGridData g) => new EdgeQuadratureScheme(true, null);
                return m_EdgeQuadraturSchemeProvider;
            }
            set {
                // deactivated due to legacy code issues:
                //if(IsCommited) 
                //    throw new NotSupportedException("not allowed to change after Commit");
                m_EdgeQuadraturSchemeProvider = value;
            }
        }

        Func<IGridData, CellQuadratureScheme> m_VolumeQuadraturSchemeProvider;
        


        /// <summary>
        /// User-customizable factory, to specify the cell/volume quadrature, see also <see cref="QuadOrderFunction"/>
        /// </summary>
        public Func<IGridData,CellQuadratureScheme> VolumeQuadraturSchemeProvider {
            get {
                if(m_VolumeQuadraturSchemeProvider == null)
                    m_VolumeQuadraturSchemeProvider = (IGridData g) => new CellQuadratureScheme(true, null);
                return m_VolumeQuadraturSchemeProvider;
            }
            set {
                //if(IsCommited)
                //    throw new NotSupportedException("not allowed to change after Commit");
                m_VolumeQuadraturSchemeProvider = value;
            }
        }

        /// <summary>
        /// Dirty hack in order to support legacy interfaces: modify quad scheme providers after commit
        /// </summary>
        internal (Func<IGridData, EdgeQuadratureScheme>,Func<IGridData,CellQuadratureScheme>) LegacySupport_ModifyQuadSchemProvider(EdgeQuadratureScheme es, CellQuadratureScheme cs) {
            var r = (EdgeQuadraturSchemeProvider, VolumeQuadraturSchemeProvider); // backup

            this.m_EdgeQuadraturSchemeProvider = g => es;
            this.m_VolumeQuadraturSchemeProvider = g => cs;

            return r;
        }

        /// <summary>
        /// Dirty hack to support legacy interface
        /// </summary>
        internal void LegacySupport_RestoreQuadSchemeProvider((Func<IGridData, EdgeQuadratureScheme> es,Func<IGridData,CellQuadratureScheme> cs) tt) {
            this.m_EdgeQuadraturSchemeProvider = tt.es;
            this.m_VolumeQuadraturSchemeProvider = tt.cs;
        }

        /// <summary>
        /// Employs <see cref="EdgeQuadraturSchemeProvider"/>, <see cref="VolumeQuadraturSchemeProvider"/>, <see cref="QuadOrderFunc"/>
        /// to generate quadrature rules for the operator evaluation.
        /// </summary>
        public (ICompositeQuadRule<QuadRule> edgeRule, ICompositeQuadRule<QuadRule> volRule) CompileQuadratureRules(IEnumerable<Basis> DomainBasis, IEnumerable<Basis> ParameterBasis, IEnumerable<Basis> CodomainBasis) {
            var order = GetOrderFromQuadOrderFunction(DomainBasis, ParameterBasis, CodomainBasis);
            IGridData gdat = DomainBasis.Any() ? DomainBasis.First().GridDat : CodomainBasis.First().GridDat;

            var edgeScheme = this.EdgeQuadraturSchemeProvider(gdat);
            var volScheme = this.VolumeQuadraturSchemeProvider(gdat);

            ICompositeQuadRule<QuadRule> _edgeRule = edgeScheme.SaveCompile(gdat, order);
            ICompositeQuadRule<QuadRule> _volRule = volScheme.SaveCompile(gdat, order);

            return (_edgeRule, _volRule);
        }

        /// <summary>
        /// internal asses for hack in <see cref="DependentTemporalOperator"/>.
        /// </summary>
        internal Dictionary<string, object> m_UserDefinedValues;

        /// <summary>
        /// Modification of <see cref="CoefficientSet.UserDefinedValues"/>, **but only if** default setting for <see cref="OperatorCoefficientsProvider"/> is used
        /// </summary>
        public IDictionary<string, object> UserDefinedValues {
            get {
                if(m_UserDefinedValues == null)
                    m_UserDefinedValues = new Dictionary<string, object>();
                return m_UserDefinedValues;
            }
            
        }


        /// <summary>
        /// <see cref="OperatorCoefficientsProvider"/>
        /// </summary>
        /// <param name="g">grid on which the operator is evaluated</param>
        /// <param name="time">current physical time</param>
        /// <returns>instance of <see cref="CoefficientSet"/> (or some derivative class)</returns>
        public delegate CoefficientSet DelOperatorCoefficientsProvider(IGridData g, double time);


        DelOperatorCoefficientsProvider m_OperatorCoefficientsProvider;
        
        
        CoefficientSet DefaultOperatorCoefficientsProvider(IGridData g, double time) {

            var r = new CoefficientSet() {
                GrdDat = g
            };

            if(g is Grid.Classic.GridData cgdat) {
                r.CellLengthScales = cgdat.Cells.CellLengthScale;
                r.EdgeLengthScales = cgdat.Edges.h_min_Edge;

            } else if(g is Grid.Aggregation.AggregationGridData agDat) { 
                r.CellLengthScales =  agDat.AncestorGrid.Cells.CellLengthScale;
                r.EdgeLengthScales =  agDat.AncestorGrid.Edges.h_min_Edge;
            } else {
                Console.Error.WriteLine("Rem: still missing cell length scales for grid type " + g.GetType().FullName);
            }

            foreach(var kv in UserDefinedValues) {
                r.UserDefinedValues[kv.Key] = kv.Value;
            }

            r.HomotopyValue = this.m_CurrentHomotopyValue;

            return r;
        }

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
                 if(IsCommitted)
                    throw new NotSupportedException("not allowed to change after operator is committed.");
                m_OperatorCoefficientsProvider = value;
            }
        }


        /// <summary>
        /// A hint for implicit/nonlinear solvers, which linearization of the operator should be used
        /// </summary>
        public LinearizationHint LinearizationHint {
            get;
            set;
        }


        static string[] GetSubarray(string[] A, int i0, int len) {
            string[] r = new string[len];
            Array.Copy(A, i0, r, 0, len);
            return r;
        }

        List<DelPartialParameterUpdate> m_ParameterUpdates = new List<DelPartialParameterUpdate>();

        /// <summary>
        /// <see cref="IDifferentialOperator.ParameterUpdates"/>
        /// </summary>
        public ICollection<DelPartialParameterUpdate> ParameterUpdates {
            get {
                if(m_IsCommitted) {
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
                if(IsCommitted) {
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
        public DifferentialOperator(int NoOfDomFields, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, params string[] __varnames)
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
        public DifferentialOperator(int NoOfDomFields, int NoOfParameters, int NoOfCodomFields, Func<int[], int[], int[], int> QuadOrderFunc, params string[] __varnames)
            : this(GetSubarray(__varnames, 0, NoOfDomFields), GetSubarray(__varnames, NoOfDomFields, NoOfParameters), GetSubarray(__varnames, NoOfDomFields + NoOfParameters, NoOfCodomFields), QuadOrderFunc) {
            if(NoOfCodomFields + NoOfDomFields + NoOfParameters != __varnames.Length)
                throw new ArgumentException("mismatch between number of provided variable names and given number of domain, parameter and codomain fields.");

        }

        /// <summary>
        /// Empty constructor; Variable, Parameter, and Codomain/Equation names are specified by the 
        /// order in which equation components are added.
        /// </summary>
        public DifferentialOperator()
            : this(new string[0], new string[0], new string[0], QuadOrderFunc.NonLinear(2)) {
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
        public DifferentialOperator(IList<string> __DomainVar, IList<string> __CoDomainVar, Func<int[], int[], int[], int> QuadOrderFunc)
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
        public DifferentialOperator(IList<string> __DomainVar, IList<string> __ParameterVar, IList<string> __CoDomainVar, Func<int[], int[], int[], int> QuadOrderFunc) {
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
            // foreach (var elem in m_ParameterVar)
            //     Console.WriteLine(elem);

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
        /// exception is thrown, if <paramref name="allowVarAddition"/> is set true, see also <see cref="IDifferentialOperator.Commit(bool)"/>
        /// </remarks>
        internal protected void Verify(bool allowVarAddition) {
            if(this.IsLinear && LinearizationHint != LinearizationHint.AdHoc)
                throw new NotSupportedException("Configuration Error: for a supposedly linear operator, the linearization hint must be " + LinearizationHint.AdHoc);

            foreach(var comps in m_EquationComponents.Values) {
                foreach(IEquationComponent c in comps) {
                    foreach(string varname in c.ArgumentOrdering) {
                        if(Array.IndexOf<string>(m_DomainVar, varname) < 0) {
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

                    if(c.ParameterOrdering != null) {
                        foreach(string varname in c.ParameterOrdering) {
                            if(Array.IndexOf<string>(m_ParameterVar, varname) < 0) {
                                if (allowVarAddition)
                                    m_ParameterVar = m_ParameterVar.Cat(varname);
                                else {
                                    throw new ApplicationException("configuration error in spatial differential operator; equation component " + c.ToString() + " depends on (parameter) variable \""
                                        + varname
                                        + "\", but this name is not a member of the parameter variable list: "
                                        + m_ParameterVar.ToConcatString("[", ", ", "]")
                                        ); ;
                                }
                            }
                        }
                    }
                }
            }


            foreach(var comps in m_EquationComponents.Values) {
                foreach(IEquationComponent c in comps) {
                    if(c.ParameterOrdering != null) {
                        foreach(string varname in c.ParameterOrdering) {
                            
                            
                            if(this.m_DomainVar.Contains(varname))
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

        _EquationComponents m_EquationComponentsHelper;

        /// <summary>
        /// for each variable in <see cref="CodomainVar"/>, a
        /// collection of equation components that define the operator.
        /// </summary>
        public IEquationComponents EquationComponents {
            get {
                return m_EquationComponentsHelper;
            }
        }

        bool m_IsCommitted = false;

        /// <summary>
        /// indicates whether the equation-assembly has been finished (by calling <see cref="Commit"/>) or not.  
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


        /// <summary>
        /// finalizes the assembly of the operator;
        /// Can be called only once in the lifetime of this object.
        /// After calling this method, no adding/removing of equation components is possible.
        /// </summary>
        public virtual void Commit(bool allowVarAddition = true) {
            Verify(allowVarAddition);

            if(TemporalOperator != null) {
                TemporalOperator.Commit();
            }

            if(m_IsCommitted)
                throw new ApplicationException("'Commit' has already been called - it can be called only once in the lifetime of this object.");

            m_IsCommitted = true;

        }


        /// <summary>
        /// implementation of <see cref="EquationComponents" />;
        /// </summary>
        public class _EquationComponents : IEquationComponents {

            internal _EquationComponents(DifferentialOperator owner) {
                m_owner = owner;
            }

            DifferentialOperator m_owner;

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
        /// returns true, if this spatial differential operator contains any 
        /// linear component,
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
        /// nonlinear component,
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

                if(Array.IndexOf<Type>(interfaces, typeof(IVolumeForm_GradUxGradV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IVolumeForm_GradUxV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IVolumeForm_UxGradV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IVolumeForm_UxV)) >= 0)
                    return true;

                if(Array.IndexOf<Type>(interfaces, typeof(IVolumeSource_GradV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IVolumeSource_V)) >= 0)
                    return true;


                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeform_GradUxGradV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeform_GradUxV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeform_UxGradV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeForm_UxV)) >= 0)
                    return true;

                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeSource_GradV)) >= 0)
                    return true;
                if(Array.IndexOf<Type>(interfaces, typeof(IEdgeSource_V)) >= 0)
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
        /// variable <paramref name="CodomVar"/> contains nonlinear integrands on edges.
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
        /// Constructs a new evaluator object for the explicit evaluation of this spatial operator.
        /// </summary>
        /// <returns></returns>
        /// <remarks>
        /// Before this method can be called,
        /// the operator assembly must be finalized by calling <see cref="Commit"/> .
        /// </remarks>
        /// <param name="CodomainVarMap">
        /// used to compute indices into the result vector
        /// </param>
        /// <param name="DomainFields">
        /// domains which are evaluated to compute fluxes, ...
        /// </param>
        /// <param name="ParameterMap">
        /// The parameter variables (of this differential operator);
        /// The number of elements in the list must match the parameter count of the differential operator
        /// (see <see cref="DifferentialOperator.ParameterVar"/>);
        /// It is allowed to set an entry to 'null', in this case the values of the parameter field
        /// are assumed to be 0.0;
        /// If the differential operator contains no parameters, this argument can be null;
        /// </param>
        public virtual IEvaluatorNonLin GetEvaluatorEx(
            IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) //
        {
            using(new FuncTrace()) {
                if(DomainFields == null)
                    DomainFields = new DGField[0];
                if(ParameterMap == null)
                    ParameterMap = new DGField[0];

                if(!IsCommitted)
                    throw new NotSupportedException("Commit() (finishing operator assembly) must be called prior to evaluation.");

                var rulz = CompileQuadratureRules(DomainFields.Select(f=>f.Basis), 
                    GetBasisS(ParameterMap),
                    CodomainVarMap.BasisS);

                var e = new EvaluatorNonLin(this, DomainFields, ParameterMap, CodomainVarMap, rulz.edgeRule, rulz.volRule);

                return e;
            }
        }

        /// <summary>
        /// Creator of a <see cref="EvaluatorLinear"/> object.
        /// </summary>
        public virtual IEvaluatorLinear GetMatrixBuilder(
            UnsetteledCoordinateMapping DomainVarMap, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) //
        {

            using(new FuncTrace()) {
                if(!IsCommitted)
                    throw new NotSupportedException("Commit() (finishing operator assembly) must be called prior to evaluation.");

                var rulz = CompileQuadratureRules((Basis[])DomainVarMap, 
                    GetBasisS(ParameterMap),
                    (Basis[])CodomainVarMap);



                var e = new EvaluatorLinear(this, DomainVarMap, ParameterMap, CodomainVarMap, rulz.edgeRule, rulz.volRule);

                return e;
            }
        }

        /// <summary>
        /// Only for debugging;  can be used to turn all edge integration in spatial operators off.
        /// </summary>
        public static bool DoEdge = true;

        /// <summary>
        /// Only for debugging; can be used to turn all volume integration in spatial operators off.
        /// </summary>
        public static bool DoVolume = true;

        /// <summary>
        /// Container for the evaluation of nonlinear fluxes/sources
        /// </summary>
        abstract public class EvaluatorBase : IEvaluator {

            DifferentialOperator m_Owner;

            /// <summary>
            /// the operator used to construct this object
            /// </summary>
            public IDifferentialOperator Owner {
                get {
                    return m_Owner;
                }
            }

            /// <summary>
            /// ctor
            /// </summary>
            protected internal EvaluatorBase(
                DifferentialOperator owner,
                UnsetteledCoordinateMapping DomainVarMap,
                IList<DGField> ParameterMap,
                UnsetteledCoordinateMapping CodomainVarMap) //
            {
                using(var tr = new FuncTrace()) {
                    MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

                    if(DomainVarMap.NoOfVariables != owner.DomainVar.Count) {
                        throw new ArgumentException("wrong number of domain variables provided.");
                    }
                    this.m_Parameters = new DGField[owner.ParameterVar.Count];
                    if(CodomainVarMap.NoOfVariables != owner.CodomainVar.Count) {
                        throw new ArgumentException("wrong number of codomain variables provided.");
                    }

                    if(!object.ReferenceEquals(DomainVarMap.GridDat, CodomainVarMap.GridDat))
                        throw new ArgumentException("Domain and Codomain map must be assigned to the same grid");

                    foreach(var f in Parameters) {
                        if(f != null) {
                            if(!object.ReferenceEquals(DomainVarMap.GridDat, f.GridDat))
                                throw new ArgumentException("Parameter fields, domain and codomain basis must be assigned to the same grid");
                        }
                    }


                    m_Owner = owner;
                    CodomainMapping = CodomainVarMap;
                    DomainMapping = DomainVarMap;
                    m_Parameters = (ParameterMap != null) ? ParameterMap.ToArray() : new DGField[0];
                    if(m_Parameters.Length != owner.ParameterVar.Count) {
                        Console.WriteLine("m_Parameters: " + m_Parameters.Length);
                        foreach (var elem in m_Parameters){
                            Console.WriteLine(elem.Identification);
                        }
                        Console.WriteLine("owner.ParameterVar: " + owner.ParameterVar.Count);
                        foreach (var elem in owner.ParameterVar){
                            Console.WriteLine(elem);
                        }
                        throw new ArgumentException("wrong number of parameter variables provided.");
                    }

                    if(!m_Owner.IsCommitted)
                        throw new ApplicationException("operator assembly must be finalized before by calling 'Commit' before this method can be called.");
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

                var _OperatorCoefficients = ((DifferentialOperator)Owner).OperatorCoefficientsProvider(this.GridData, this.time);


                Debug.Assert(CodNames.Length == CodDGdeg.Length);
                for(int iCod = 0; iCod < CodDGdeg.Length; iCod++) {
                    var comps = m_Owner.m_EquationComponents[CodNames[iCod]];
                    foreach(var c in comps) {
                        if(c is IEquationComponentCoefficient) {
                            var ce = c as IEquationComponentCoefficient;

                            int[] DomDGdeg_cd = new int[ce.ArgumentOrdering.Count];
                            for(int i = 0; i < DomDGdeg_cd.Length; i++) {
                                string domName = ce.ArgumentOrdering[i];
                                int idx = Array.IndexOf(DomNames, domName);
                                DomDGdeg_cd[i] = DomDGdeg[idx];
                            }


                            ce.CoefficientUpdate(_OperatorCoefficients, DomDGdeg_cd, CodDGdeg[iCod]);
                        }
                    }
                }
            }


            SubGridBoundaryModes m_SubGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge;

            /// <summary>
            /// 
            /// </summary>
            protected CellMask m_SubGrid_InCells;


            /// <summary>
            /// State set by <see cref="ActivateSubgridBoundary"/>
            /// </summary>
            public SubGridBoundaryModes SubGridBoundaryTreatment {
                get {
                    return m_SubGridBoundaryTreatment;
                }
            }

            /// <summary>
            /// Restricts the evaluation of the operator to a specific cell mask.
            /// </summary>
            /// <param name="Mask">
            /// cell mask where the operator should be evaluated
            /// </param>
            /// <param name="subGridBoundaryTreatment">
            /// defines what is to be done at edges where one neighboring cell is part of the cell mask <paramref name="Mask"/>, 
            /// but the other neighboring cell is *not* part of the cell mask <paramref name="Mask"/>.
            /// </param>
            public void ActivateSubgridBoundary(CellMask Mask, SubGridBoundaryModes subGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge) {
                if(!object.ReferenceEquals(Mask.GridData, this.GridData))
                    throw new ArgumentException("grid mismatch");
                if(Mask != null && Mask.MaskType != MaskType.Logical)
                    throw new ArgumentException("expecting logical mask");
                m_SubGrid_InCells = Mask;
                m_SubGridBoundaryTreatment = subGridBoundaryTreatment;
            }


            /// <summary>
            /// coordinate mapping for the codomain variables;
            /// </summary>
            public UnsetteledCoordinateMapping CodomainMapping {
                get;
                private set;
            }



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
                get;
                private set;
            }

            /// <summary>
            /// Grid, on which this evaluator operates on.
            /// </summary>
            public IGridData GridData {
                get {
                    Debug.Assert(object.ReferenceEquals(DomainMapping.GridDat, CodomainMapping.GridDat));
                    return DomainMapping.GridDat;
                }
            }




            double m_time = 0.0;

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


        internal bool RequiresComplicatedPeriodicity(IGridData gd) {
            /*
            var tags = gd.iGeomEdges.EdgeTags;
            int E = gd.iGeomEdges.Count;
            Debug.Assert(E == tags.Length);
            
            //for(int e = 0; e < E; e++) {
            //    int tag = tags[e];
            //    if(tag >= Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG) {

            //        return true;
            //    }
            //}

            foreach(var Trafo in gd.Grid.PeriodicTrafo) {
                var Mtx = Trafo.Matrix.CloneAs();
                Mtx.AccEye(-1.0);
                if(Mtx.InfNorm() >= 1e-8)
                    throw new NotSupportedException("Non-parallel periodic edges are still not supported.");
            }
            */

            // feature deactivated for now, because the implementation for linear components sucks so much

            return false;
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
            /// If the grid contains some periodic boundaries which are not parallel (e.g. some cake-pie-subsection of a rotational domain)
            /// periodicity required additional transformations/rotations of vectors for **both** sides of the periodic edge;
            /// Furthermore, these rotations are different (inverse) for the in- and the out-edge, 
            /// therefore the contribution to the out-cell is computed in a second pass, by this integrator.
            /// </summary>
            BoSSS.Foundation.Quadrature.NonLin.NECQuadratureEdge m_ComplicatedPeriodicEdge;

            /// <summary>
            /// if the right-hand-side is present and contains nonlinear components, 
            /// this is the corresponding volume quadrature;
            /// otherwise, this member is null;
            /// </summary>
            BoSSS.Foundation.Quadrature.NonLin.NECQuadratureVolume m_NonlinearVolume;

            // helper hack to support empty domain lists
            static UnsetteledCoordinateMapping Helper(
                IList<DGField> DomainVarMap,
                IList<DGField> ParameterMap) {
                IGridData g = null;
                if(DomainVarMap != null && DomainVarMap.Count > 0) {
                    g = DomainVarMap[0].Basis.GridDat;
                } else if(ParameterMap != null && ParameterMap.Count > 0) {
                    g = ParameterMap[0].Basis.GridDat;
                }

                if(DomainVarMap != null && DomainVarMap.Count > 0) {
                    return new UnsetteledCoordinateMapping(DomainVarMap.Select(f => f.Basis));
                } else {
                    return new UnsetteledCoordinateMapping(g);
                }
            }


            /// <summary>
            /// Not for direct user interaction
            /// </summary>
            internal protected EvaluatorNonLin(
                DifferentialOperator owner,
                IList<DGField> DomainVarMap,
                IList<DGField> ParameterMap,
                UnsetteledCoordinateMapping CodomainVarMap,
                ICompositeQuadRule<QuadRule> edgeQuadRule,
                ICompositeQuadRule<QuadRule> volQuadRule) //
             : base(owner, Helper(DomainVarMap, ParameterMap), ParameterMap, CodomainVarMap) //
            {

                var grdDat = base.GridData;
                if(DomainVarMap != null & DomainVarMap.Count > 0)
                    DomainFields = new CoordinateMapping(DomainVarMap);
                else
                    DomainFields = new CoordinateMapping(grdDat);

                

                if(owner.RequiresEdgeQuadrature) {


                    m_NonlinearEdge = new BoSSS.Foundation.Quadrature.NonLin.NECQuadratureEdge(grdDat,
                                                            (DifferentialOperator) Owner,
                                                            DomainVarMap,
                                                            ParameterMap,
                                                            CodomainVarMap,
                                                            edgeQuadRule);

                    if(owner.RequiresComplicatedPeriodicity(grdDat)) {
                        m_ComplicatedPeriodicEdge = new BoSSS.Foundation.Quadrature.NonLin.NECQuadratureEdge(grdDat,
                                                            (DifferentialOperator)Owner,
                                                            DomainVarMap,
                                                            ParameterMap,
                                                            CodomainVarMap,
                                                            RestrictQr(edgeQuadRule, GetPeriodicEdgesMask(grdDat)));

                        m_NonlinearEdge._PeriodicVectorTrafo = Foundation.Quadrature.NonLin.NECQuadratureEdge.PeriodicVectorTrafo.bck;
                        m_ComplicatedPeriodicEdge._PeriodicVectorTrafo = Foundation.Quadrature.NonLin.NECQuadratureEdge.PeriodicVectorTrafo.fwd;
                    }


                }

                if(owner.RequiresVolumeQuadrature) {
                    m_NonlinearVolume = new BoSSS.Foundation.Quadrature.NonLin.NECQuadratureVolume(grdDat,
                                                                (DifferentialOperator) Owner,
                                                                DomainVarMap,
                                                                ParameterMap,
                                                                CodomainVarMap,
                                                                volQuadRule);


                }


                base.MPITtransceive = true;
            }

            static EdgeMask GetPeriodicEdgesMask(IGridData gd) {
                var tags = gd.iGeomEdges.EdgeTags;
                int E = gd.iGeomEdges.Count;
                Debug.Assert(E == tags.Length);
                var bmsk = new BitArray(E);
                for(int e = 0; e < E; e++) {
                    if(tags[e] >= Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG)
                        bmsk[e] = true;
                }
                return new EdgeMask(gd, bmsk);
            }

            /// <summary>
            /// Removes all edes from a quadrature rule <paramref name="Qr"/> which are not contained in the mask <paramref name="Restriction"/>
            /// </summary>
            /// <remarks>
            /// Brute-force implementation, not very sophisticated; should be ok for now.
            /// </remarks>
            static ICompositeQuadRule<QuadRule> RestrictQr(ICompositeQuadRule<QuadRule> Qr, EdgeMask Restriction) {
                
                int E = Restriction.GridData.iGeomEdges.Count;
                var RestrictionBMask = Restriction.GetBitMask();
                Debug.Assert(RestrictionBMask.Length == E);

                QuadRule[] temp = new QuadRule[E];
                foreach(var pair in Qr) {
                    for(int i = pair.Chunk.i0; i < pair.Chunk.JE; i++) {
                        if(RestrictionBMask[i])
                            temp[i] = pair.Rule;
                    }
                }

                var fQr = new CompositeQuadRule<QuadRule>();
                for(int i = 0; i < E; i++) {
                    if(temp[i] != null) {
                        int iLast = fQr.chunkRulePairs.Count - 1;
                        if(fQr.chunkRulePairs.Count > 0
                            && fQr.chunkRulePairs[iLast].Chunk.JE == i
                            && object.ReferenceEquals(fQr.chunkRulePairs[iLast].Rule, temp[i])) {
                            var lastShit = fQr.chunkRulePairs[iLast].Chunk;
                            var enlargedChunk = new Chunk();
                            enlargedChunk.i0 = lastShit.i0;
                            enlargedChunk.Len = lastShit.Len + 1;

                            var newShit = new ChunkRulePair<QuadRule>(enlargedChunk, temp[i]);
                            fQr.chunkRulePairs[iLast] = newShit;
                        } else {
                            var next = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(i), temp[i]);
                            fQr.chunkRulePairs.Add(next);
                        }

                        Debug.Assert(fQr.chunkRulePairs.Last().Chunk.JE == i + 1);
                    }

                }

#if DEBUG
                {
                    QuadRule[] check = new QuadRule[E];
                    foreach(var pair in fQr) {
                        for(int i = pair.Chunk.i0; i < pair.Chunk.JE; i++) {
                            if(RestrictionBMask[i])
                                check[i] = pair.Rule;
                        }
                    }

                    for(int i = 0; i < E; i++) {
                        if(check[i] != null) {
                            Debug.Assert(RestrictionBMask[i] == true);
                            Debug.Assert(object.ReferenceEquals(check[i], temp[i]));
                        } else {
                            Debug.Assert(RestrictionBMask[i] == false || temp[i] == null);
                        }
                    }

                }

#endif

                return fQr;
            }



            /// <summary>
            /// DG fields which serve a input for the spatial operator.
            /// </summary>
            public CoordinateMapping DomainFields {
                get;
                private set;
            }

            /// <summary>
            /// evaluates the differential operator (<see cref="EvaluatorBase.Owner"/>)
            /// for the domain variables/fields in <see cref="EvaluatorBase.DomainMapping"/>, i.e.
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
            /// Indices into this vector are computed according to <see cref="EvaluatorBase.CodomainMapping"/>;
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

                    if(m_NonlinearVolume != null && DoVolume) {
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



                    void CallEdge(Quadrature.NonLin.NECQuadratureEdge ne, string name) {
                        if(ne != null && DoEdge) {
                            using(new BlockTrace(name, tr)) {

                                ne.m_Output = output;
                                ne.m_alpha = alpha;
                                ne.Time = time;
                                ne.SubGridBoundaryTreatment = base.SubGridBoundaryTreatment;
                                ne.SubGridCellsMarker = (base.m_SubGrid_InCells != null) ? base.m_SubGrid_InCells.GetBitMaskWithExternal() : null;

                                ne.m_outputBndEdge = outputBndEdge;

                                ne.Execute();

                                ne.m_Output = null;
                                ne.m_outputBndEdge = null;
                                ne.m_alpha = 1.0;
                                ne.SubGridCellsMarker = null;
                            }
                        }
                    }

                    CallEdge(m_NonlinearEdge, "Edge_Integration_NonLin");
                    CallEdge(m_ComplicatedPeriodicEdge, "Edge_Integration_NonLin_periodic");


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
                DifferentialOperator owner,
                UnsetteledCoordinateMapping DomainVarMap,
                IList<DGField> ParameterMap,
                UnsetteledCoordinateMapping CodomainVarMap,
                ICompositeQuadRule<QuadRule> edgeQr,
                ICompositeQuadRule<QuadRule> volQr) //
                 : base(owner, DomainVarMap, ParameterMap, CodomainVarMap) //
            {
                owner.RequiresComplicatedPeriodicity(CodomainMapping.GridDat);
                
                foreach(string codVarName in owner.CodomainVar) {
                    var comps = owner.EquationComponents[codVarName];

                    //if (comps.Where(cmp => cmp is INonlinearFlux).Count() > 0)
                    //    throw new NotSupportedException("'INonlinearFlux' is not supported for linearization; (codomain variable '" + codVarName + "')");
                    if(comps.Where(cmp => cmp is INonlinearFluxEx).Count() > 0)
                        throw new NotSupportedException("'INonlinearFluxEx' is not supported for linearization; (codomain variable '" + codVarName + "')");
                    //if (comps.Where(cmp => cmp is IDualValueFlux).Count() > 0)
                    //    throw new NotSupportedException("'IDualValueFlux' is not supported for linearization; (codomain variable '" + codVarName + "')");
                    if(comps.Where(cmp => cmp is INonlinearSource).Count() > 0)
                        throw new NotSupportedException("'INonlinearSource' is not supported for linearization; (codomain variable '" + codVarName + "')");
                }

                this.edgeRule = edgeQr;
                this.volRule = volQr;
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
                Internal_ComputeMatrixEx(default(BlockMsrMatrix), AffineOffset, true, 1.0);
            }

            /// <summary>
            /// computes a linearization of the operator in the form 
            /// \f[
            ///    \mathcal{M} U + \mathcal{B}.
            /// \f]
            /// </summary>
            /// <param name="Matrix">
            /// Output, the operator matrix, scaled by <paramref name="alpha"/>, is accumulated here
            /// </param>
            /// <param name="AffineOffset">
            /// Output, the affine part of the operator linearization, scaled by <paramref name="alpha"/>, is accumulated here
            /// </param>
            /// <param name="alpha">
            /// scaling factor
            /// </param>
            public void ComputeMatrix<M, V>(M Matrix, V AffineOffset, double alpha = 1.0)
                where M : IMutableMatrixEx
                where V : IList<double> // 
            {
                Internal_ComputeMatrixEx(Matrix, AffineOffset, false, alpha);
            }

            /// <summary>
            /// returns parameter fields
            /// </summary>
            protected override DGField[] GetTrxFields() {
                if(Parameters == null) {
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
                M Matrix, V AffineOffset, bool OnlyAffine, double alpha)
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
                    

                    
                    DifferentialOperator _Owner = (DifferentialOperator)this.Owner;
                    
                    if(volRule.Any() && DoVolume) {
                        using(new BlockTrace("Volume_Integration_(new)", tr)) {

                            /*
                            if (Matrix != null && AffineOffset != null) {
                                for (int i = 0; i < 1; i++) {
                                    Console.WriteLine($"   i = {i}");

                                    double afferrNorm = 0, mtxerrNorm = 0;
                                    for (int k = 0; k < 1; k++) {
                                        var clMatrix = new MsrMatrix(Matrix.RowPartitioning, Matrix.ColPartition);
                                        var clAffine = new double[AffineOffset.Count];
                                        var mpMatrix = new MsrMatrix(Matrix.RowPartitioning, Matrix.ColPartition);
                                        var mpAffine = new double[AffineOffset.Count];

                                        Debugi.SkipComp = -1;
                                        Debugi.CompCont = 0;
                                        Debugi.printInfo = k == 0;
                                        ilPSP.Environment.NumThreads = 1;
                                        var clmtxBuilder = new LECVolumeQuadrature2<MsrMatrix, double[]>(_Owner);
                                        clmtxBuilder.m_alpha = alpha;
                                        clmtxBuilder.Execute(volRule, CodomainMapping, Parameters, DomainMapping, clMatrix, clAffine, time);

                                        Debugi.SkipComp = -1;
                                        Debugi.CompCont = 0;
                                        //Debugi.printInfo = false;
                                        ilPSP.Environment.NumThreads = 8;
                                        var mpmtxBuilder = new LECVolumeQuadrature2<MsrMatrix, double[]>(_Owner);
                                        mpmtxBuilder.m_alpha = alpha;
                                        mpmtxBuilder.Execute(volRule, CodomainMapping, Parameters, DomainMapping, mpMatrix, mpAffine, time);


                                        var errMtx = clMatrix.CloneAs();
                                        errMtx.Acc(mpMatrix, -1.0);
                                        mtxerrNorm += errMtx.InfNorm();

                                        afferrNorm += clAffine.MPI_L2Dist(mpAffine);

                                    }
                                    Console.WriteLine($"   difference (i = {i}): {mtxerrNorm:0.####e-00}  {mtxerrNorm:0.####e-00}");
                                    Console.Write("");
                                }
                            }
                            //*/

                            var mtxBuilder = new LECVolumeQuadrature2<M, V>(_Owner);
                            mtxBuilder.m_alpha = alpha;
                            mtxBuilder.Execute(volRule, CodomainMapping, Parameters, DomainMapping, OnlyAffine ? default(M) : Matrix, AffineOffset, time);

                            //volRule.ToTextFileVolume(this.GridData as BoSSS.Foundation.Grid.Classic.GridData, "Volume.csv");
                        }

                    } else {
                        //tr.Info("volume integration skipped: cell mask is empty");
                    }
                    //*/

                    // edge integration
                    // ----------------
                    
                    if(!edgeRule.IsNullOrEmpty() && DoEdge) {
                        using(new BlockTrace("Edge_Integration_(new)", tr)) {
                            var mxtbuilder2 = new LECEdgeQuadrature2<M, V>(_Owner);
                            mxtbuilder2.m_alpha = alpha;
                            mxtbuilder2.Execute(edgeRule, CodomainMapping, Parameters, DomainMapping, OnlyAffine ? default(M) : Matrix, AffineOffset, time);
                        }
                    }
                }
            }
        }


        /// <summary>
        /// constructs a <see cref="FDJacobianBuilder"/> object to linearize nonlinear operators
        /// </summary>
        public virtual IEvaluatorLinear GetFDJacobianBuilder(
            IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap) //
        {

            Action<double, IEnumerable<DGField>, IEnumerable<DGField>> ParamUpdate =
                delegate (double time, IEnumerable<DGField> DomF, IEnumerable<DGField> ParamF) {
                    this.InvokeParameterUpdate(time, DomF.ToArray(), ParamF.ToArray());
                };

            return GetFDJacobianBuilder_(DomainFields, ParameterMap, CodomainVarMap, ParamUpdate);
        }

        static Basis[] GetBasisS(IList<DGField> ParameterMap) {
            if(ParameterMap == null)
                return new Basis[0];

            return ParameterMap.Select(f => f != null ? f.Basis : default(Basis)).ToArray();
        }

        /// <summary>
        /// Internal implementation and legacy API;
        /// constructs a <see cref="FDJacobianBuilder"/> object to linearize nonlinear operators
        /// </summary>
        /// <param name="CodomainVarMap"></param>
        /// <param name="DomainFields"></param>
        /// <param name="ParameterMap"></param>
        /// <param name="legayc_delParameterUpdate">
        /// legacy: external delegate to update all parameters at once;
        /// specifying this replaces all <see cref="IParameterHandling"/> components and all <see cref="ParameterUpdates"/> set for this operator
        /// with the external update.
        /// </param>
        public virtual FDJacobianBuilder GetFDJacobianBuilder_(
            IList<DGField> DomainFields, IList<DGField> ParameterMap, UnsetteledCoordinateMapping CodomainVarMap,
            Action<double, IEnumerable<DGField>, IEnumerable<DGField>> legayc_delParameterUpdate) //
        {
            using(new FuncTrace()) {
                if(!IsCommitted)
                    throw new NotSupportedException("Commit() (finishing operator assembly) must be called prior to evaluation.");

                
                var rulz = CompileQuadratureRules(DomainFields.Select(f=>f.Basis), 
                    GetBasisS(ParameterMap),
                    CodomainVarMap.BasisS);


                var e = new FDJacobianBuilder(new EvaluatorNonLin(
                    this,
                    new CoordinateMapping(DomainFields), ParameterMap, CodomainVarMap,
                    rulz.edgeRule, rulz.volRule),
                    legayc_delParameterUpdate);
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
            public FDJacobianBuilder(IEvaluatorNonLin __Eval, Action<double, IEnumerable<DGField>, IEnumerable<DGField>> __delParameterUpdate) {

                eps = 1.0;
                while(1.0 + eps > 1.0) {
                    eps = eps / 2;
                }
                eps = Math.Sqrt(eps);

                Eval = __Eval;
                Eval.MPITtransceive = true;
                if(__delParameterUpdate != null)
                    DelParamUpdate = __delParameterUpdate;
                else
                    DelParamUpdate = EmptyParamUpdate;

                BuildOptimizedGridColoring();
                //BuildOneByOneColoring();

                //Console.WriteLine("FDJac: no of color lists: " + ColorLists.Length);
            }

            void EmptyParamUpdate(double t, IEnumerable<DGField> F, IEnumerable<DGField> P) {
                // do nothing
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
            public IDifferentialOperator Owner {
                get {
                    return Eval.Owner;
                }
            }

            /// <summary>
            /// can be dangerous to turn off
            /// </summary>
            public bool MPITtransceive {
                get {
                    return Eval.MPITtransceive;
                }
                set {
                    Eval.MPITtransceive = value;
                }
            }

            /// <summary>
            /// Internally used evaluation for finite differences
            /// </summary>
            public IEvaluatorNonLin Eval {
                get;
                private set;
            }

            Action<double, IEnumerable<DGField>, IEnumerable<DGField>> DelParamUpdate;

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

            /// <summary>
            /// coloring which exploits presumed locality of the DG operator
            /// </summary>
            void BuildOptimizedGridColoring() {
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
                        if(Colored[j] == true)
                            continue;
                        if(LocalMarker[j] != 0)
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

                            for(int je = J; je < JE; je++) {
                                if(LocalMarker[je] != 0 && ExchangedMarker[je] != 0) {
                                    Debug.Assert(LocalMarker[je] != ExchangedMarker[je]);
                                    LocalConflicts++;

                                    double rndVal = rnd.NextDouble();

                                    //// some parallel conflict detected: one of the two ranks has to yield

                                    //if (ExchangedMarker[je] > 0 && Math.Abs(ExchangedMarker[je]) > myMarkerToken) {
                                    //    // the other rank should yield
                                    //} else {
                                    //    // this rank has to yield
                                    if(rndVal >= 0.5) {
                                        int jToRemove = LocalColorCause[je];
                                        Debug.Assert(jToRemove < J);
                                        Debug.Assert(ColoredPass[jToRemove] == true);

                                        Removed.Add(jToRemove);

                                        ColoredPass[jToRemove] = false;
                                        LocalMarker[jToRemove] = 0;
                                        int[] Neighs_jToRemove = Neighs[jToRemove];
                                        foreach(int jn in Neighs_jToRemove) {
                                            LocalMarker[jn] = 0;
                                        }

                                        //}
                                    }
                                }
                            }

                            GlobalConflicts = LocalConflicts.MPISum();

                        } while(GlobalConflicts > 0);

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
                    if(GlobColoredPass <= 0) {
                        DeadlockWatch++;
                        if(DeadlockWatch >= 1000)
                            throw new ApplicationException("Deadlock in parallel coloring.");
                        continue;
                        // dbg_launch();
                    }


                    ColorListsTmp.Add(CellList.ToArray());

                    // communicate external lists
                    // ==========================

                    if(gDat.MpiSize > 1) {

                        // dbg_launch();

                        var ExchData = new Dictionary<int, List<Tuple<long, long>>>();

                        foreach(int j in CellList) {
                            int[] Neighs_j = Neighs[j];
                            foreach(int jN in Neighs_j) {
                                if(jN >= J) {

                                    long Gl_jN = GlidxExt[jN - J];
                                    int iProc = CellPart.FindProcess(Gl_jN);
                                    long Gl_j = j + CellPart.i0;

                                    if(!ExchData.TryGetValue(iProc, out var ExchData_iProc)) {
                                        ExchData_iProc = new List<Tuple<long, long>>();
                                        ExchData.Add(iProc, ExchData_iProc);
                                    }

                                    ExchData_iProc.Add(new Tuple<long, long>(Gl_j, Gl_jN));
                                }
                            }
                        }

                        var RcvData = SerialisationMessenger.ExchangeData(ExchData);

                        var ExtColor = new Dictionary<int, List<int>>();

                        foreach(var kv in RcvData) {
                            int iProc = kv.Key;
                            var list = kv.Value;

                            foreach(var t in list) {
                                long Gl_j = t.Item1;
                                long Gl_jN = t.Item2;
                                Debug.Assert(CellPart.FindProcess(Gl_j) == iProc);
                                Debug.Assert(CellPart.IsInLocalRange(Gl_jN));

                                int Loc_jN = checked((int)(Gl_jN - CellPart.i0));
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
            /// coloring with one cell at a pass
            /// </summary>
            void BuildOneByOneColoring() {
                var gDat = Eval.GridData;

                int[][] Neighs = gDat.iLogicalCells.CellNeighbours;
                int J = gDat.iLogicalCells.NoOfLocalUpdatedCells;
                int JE = gDat.iLogicalCells.Count;
                long[] GlidxExt = gDat.iParallel.GlobalIndicesExternalCells;
                var Gl2LocExt = gDat.iParallel.Global2LocalIdx;
                var CellPart = gDat.CellPartitioning;
                long Jglob = CellPart.TotalLength;


                //int[] LocalMarker = new int[JE]; //    marker for blocked in the current pass 
                //int[] ExchangedMarker = new int[JE]; //  accumulation buffer for MPI exchange
                //BitArray Colored = new BitArray(JE); // all cells which are already colored (in previous passes)
                //BitArray ColoredPass = new BitArray(JE); // all cells which are colored in current pass
                //int[] LocalColorCause = new int[JE];

                List<int> CellList = new List<int>();
                List<int[]> ColorListsTmp = new List<int[]>();
                List<int[]> ExternalColorListsTmp = new List<int[]>();
                List<int[][]> ExternalColorListsNeighborsTmp = new List<int[][]>();

                for(int jClGlob = 0; jClGlob < Jglob; jClGlob++) {

                    // find next color list
                    // ====================

                    int jLoc = -1234;
                    if(CellPart.IsInLocalRange(jClGlob)) {
                        jLoc = CellPart.TransformIndexToLocal(jClGlob);
                        ColorListsTmp.Add(new int[] { jLoc });
                    }


                    // communicate external lists
                    // ==========================

                    if(gDat.MpiSize > 1) {

                        // dbg_launch();

                        var ExchData = new Dictionary<int, List<Tuple<int, int>>>();

                        if(jLoc >= 0) {
                            int[] Neighs_j = Neighs[jLoc];
                            foreach(int jN in Neighs_j) {
                                if(jN >= J) {

                                    int Gl_jN = (int)GlidxExt[jN - J];
                                    int iProc = CellPart.FindProcess(Gl_jN);
                                    int Gl_j = jClGlob;

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

                        foreach(var kv in RcvData) {
                            int iProc = kv.Key;
                            var list = kv.Value;

                            foreach(var t in list) {
                                int Gl_j = t.Item1;
                                int Gl_jN = t.Item2;
                                Debug.Assert(CellPart.FindProcess(Gl_j) == iProc);
                                Debug.Assert(CellPart.IsInLocalRange(Gl_jN));

                                int Loc_jN = checked((int)(Gl_jN - CellPart.i0));
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
            /// Output, the approximate Jacobian matrix of the operator, scaled by <paramref name="alpha"/>, is accumulated here
            /// </param>
            /// <param name="AffineOffset">
            /// Output, the operator value in the linearization point, scaled by <paramref name="alpha"/>.
            /// </param>
            /// <param name="alpha">
            /// scaling factor
            /// </param>
            public void ComputeMatrix<M, V>(M Matrix, V AffineOffset, double alpha = 1.0)
                where M : IMutableMatrixEx
                where V : IList<double> // 
            {

                // init locals
                // ===========
                var codMap = Eval.CodomainMapping;
                var domMap = Eval.DomainMapping;
                DGField[] domFields = Eval.DomainFields.Fields.ToArray();
                var U0 = new CoordinateVector(Eval.DomainFields);

                long j0 = Eval.GridData.CellPartitioning.i0;
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
                DelParamUpdate(this.time, domFields, Eval.Parameters.ToArray());
#if DEBUG
                CoordinateVector ParamsVec;
                double[] ParamsVecBkup;
                if (Eval.Parameters.Where(f => f != null).Count() > 0) {
                    var AllocatedParams = Eval.Parameters.Where(f => f != null).ToArray();
                    ParamsVec = new CoordinateVector(AllocatedParams);
                    ParamsVecBkup = ParamsVec.ToArray();
                } else {
                    ParamsVec = null;
                    ParamsVecBkup = null;
                }
#endif
                Eval.Evaluate(1.0, 0.0, F0);
                AffineOffset.AccV(alpha, F0);
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

                for(int iCellPass = 0; iCellPass < ColorLists.Length; iCellPass++) { // loop over all cell lists...
                    int[] CellList = this.ColorLists[iCellPass];
                    int[] ExtCellList = this.ExternalColorLists[iCellPass];
                    
                    int[] CoordCounter = new int[JE];
                    int[] FieldCounter = new int[JE];

                    int maxNj = 0;
                    foreach(int j in CellList) {
                        int Nj = domMap.GetTotalNoOfCoordinatesPerCell(j);
                        maxNj = Math.Max(Nj, maxNj);
                    }
                    maxNj = maxNj.MPIMax();

                    Buffer.Clear();

                    for(int n = 0; n < maxNj; n++) { // loop over DG coordinates in cell

                        // backup DG coordinates
                        // ---------------------
                        U0backup.SetV(U0);

                        // apply distortions
                        // -----------------
                        int AnyLoc = 0;
                        foreach(int j in CellList) {
                            int iFld = FieldCounter[j];
                            int nFld = CoordCounter[j];
                            if(iFld >= NoOfDomFields)
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
                        if(AnyGlob >= 0)
                            break; // finished with entire cell list on all processors

                        // evaluate operator
                        // -------------------
                        EvalBuf.ClearEntries();
                        DelParamUpdate(this.time, domFields, Eval.Parameters.ToArray());
                        Eval.Evaluate(1.0, 0.0, EvalBuf);
                        NoOfEvals++;

                        // ------------------------------

                        for(int IntExt = 0; IntExt < 2; IntExt++) {
                            int[] __CellList;
                            switch(IntExt) {
                                case 0: __CellList = CellList; break;
                                case 1: __CellList = ExtCellList; break;
                                default: throw new ApplicationException();
                            }


                            // save results
                            // -------------------------------
                            int cnt = 0;
                            foreach(int _j in __CellList) {
                                int[] Neighs_j; // = Neighs[_j];
                                switch(IntExt) {
                                    case 0: Neighs_j = Neighs[_j]; break;
                                    case 1: Neighs_j = this.ExternalColorListsNeighbors[iCellPass][cnt]; break;
                                    default: throw new ApplicationException();
                                }
                                cnt++;

                                int jCol = _j;

                                int iFldCol = FieldCounter[jCol];
                                int nFldCol = CoordCounter[jCol];
                                if(iFldCol >= NoOfDomFields)
                                    continue; // finished with cell

                                int iCol = domMap.LocalUniqueCoordinateIndex(iFldCol, jCol, nFldCol);
                                int i0Col = domMap.LocalUniqueCoordinateIndex(0, jCol, 0);
                                int iRelCol = iCol - i0Col;

                                for(int k = 0; k <= Neighs_j.Length; k++) { // loop over neighbors which are influenced by the distortion
                                    int jRow;
                                    if(k == 0) {
                                        jRow = _j;
                                    } else {
                                        jRow = Neighs_j[k - 1];
                                    }

                                    if(jRow >= J) {
                                        continue; // external cell; should be treated on other proc.
                                    }
                                    int i0Row = codMap.LocalUniqueCoordinateIndex(0, jRow, 0);
                                    int NoOfRows = codMap.GetBlockLen(jRow);

                                    for(int iRelRow = 0; iRelRow < NoOfRows; iRelRow++) {
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
                            foreach(int j in __CellList) {
                                int iFld = FieldCounter[j];
                                if(iFld >= NoOfDomFields)
                                    continue; // finished with cell 'j'

                                int Nj = domMap.BasisS[iFld].GetLength(j);
                                CoordCounter[j]++;
                                if(CoordCounter[j] >= Nj) {
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

                    for(int IntExt = 0; IntExt < 2; IntExt++) {
                        int[] __CellList;
                        switch(IntExt) {
                            case 0: __CellList = CellList; break;
                            case 1: __CellList = ExtCellList; break;
                            default: throw new ApplicationException();
                        }

                        int cnt = 0;
                        foreach(int _j in __CellList) {
                            int[] Neighs_j; // = Neighs[_j];
                            switch(IntExt) {
                                case 0: Neighs_j = Neighs[_j]; break;
                                case 1: Neighs_j = this.ExternalColorListsNeighbors[iCellPass][cnt]; break;
                                default: throw new ApplicationException();
                            }
                            cnt++;

                            int jCol = _j;
                            int i0Col = domMap.LocalUniqueCoordinateIndex(0, jCol, 0);
                            int iECol = domMap.LocalUniqueCoordinateIndex(NoOfDomFields - 1, jCol, lastDomB.GetLength(jCol) - 1);

                            for(int k = 0; k <= Neighs_j.Length; k++) { // loop over neighbors which are influenced by the distortion
                                int jRow;
                                if(k == 0) {
                                    jRow = _j;
                                } else {
                                    jRow = Neighs_j[k - 1];
                                }

                                if(jRow >= J)
                                    continue; // external cell; should be treated on other proc.


                                int i0Row = codMap.LocalUniqueCoordinateIndex(0, jRow, 0);
                                int iERow = codMap.LocalUniqueCoordinateIndex(NoOfCodFields - 1, jRow, lastCodB.GetLength(jRow) - 1);

                                var Block = Buffer.ExtractSubArrayShallow(new int[] { i0Row, 0 }, new int[] { iERow, iECol - i0Col });

                                Matrix.AccBlock(i0Row + codMap.i0,
                                    //i0Col + domMap.i0, 
                                    domMap.GlobalUniqueCoordinateIndex(0, jCol, 0),
                                    alpha, Block);
                            }
                        }

                    }
                }
                // restore original state before return
                // ====================================
                U0.SetV(U0backup);
                DelParamUpdate(this.time, domFields, Eval.Parameters.ToArray());
#if DEBUG
                if (Eval.Parameters.Count > 0) {
                    double deltaParamsVec = ParamsVecBkup.L2DistPow2(ParamsVecBkup).MPISum().Sqrt();
                    double abs = Math.Max(Math.Sqrt(BLAS.MachineEps), Math.Max(ParamsVecBkup.MPI_L2Norm(), ParamsVec.MPI_L2Norm()));
                    double relDelta = deltaParamsVec / abs;
                    if (deltaParamsVec > Math.Sqrt(BLAS.MachineEps))
                        throw new ApplicationException("FDJacobian detected side-effect in DelParamUpdate(...) -- provides different result before and after finite difference approximation.");
                }

#endif

                //Console.WriteLine("Total number of evaluations: " + NoOfEvals);

                // correct the affine offset
                // =========================

                // actually, the Jacobian is the approx M*(U-U0) + b
                // but we want                          M*U + (b - M*U0)
                // (in this fashion, this approximation also works with BDF schemes *in the same fashion* as other linearizations
                Matrix.SpMV(-1.0, U0, 1.0, AffineOffset);
            }

            /// <summary>
            /// Evaluation at the linearization point
            /// </summary>
            public void ComputeAffine<V>(V AffineOffset) where V : IList<double> {
                //int Lout = Eval.CodomainMapping.LocalLength;

                //double[] F0 = new double[Lout];
                //DelParamUpdate(Eval.DomainFields.Fields.ToArray(), Eval.Parameters.ToArray());
                //Eval.Evaluate(1.0, 0.0, F0);
                //AffineOffset.AccV(1.0, F0);

                BlockMsrMatrix dummy = new BlockMsrMatrix(this.CodomainMapping, this.DomainMapping);
                this.ComputeMatrix(dummy, AffineOffset);
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="sgrd"></param>
            /// <param name="subGridBoundaryTreatment"></param>
            public void ActivateSubgridBoundary(CellMask sgrd, SubGridBoundaryModes subGridBoundaryTreatment = SubGridBoundaryModes.BoundaryEdge) {
                if(sgrd != null && sgrd.MaskType != MaskType.Logical)
                    throw new ArgumentException("expecting logical mask");
                Eval.ActivateSubgridBoundary(sgrd, subGridBoundaryTreatment);
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
        public DifferentialOperator _GetJacobiOperator(int SpatialDimension) {
            if(!this.IsCommitted)
                throw new InvalidOperationException("Invalid prior to calling Commit().");

            // parameters and activation flags
            // ===============================

            var allcomps = new List<IEquationComponent>();
            foreach(var cdo in this.CodomainVar)
                allcomps.AddRange(this.EquationComponents[cdo]);

            TermActivationFlags extractTaf(IEquationComponent c) {
                TermActivationFlags ret = default(TermActivationFlags);
                if(c is IVolumeForm vf) {
                    ret = ret | vf.VolTerms;
                }

                if(c is IEdgeForm ef) {
                    ret = ret | ef.BoundaryEdgeTerms;
                    ret = ret | ef.InnerEdgeTerms;
                }

                return ret;
            }

            var h = new JacobianParamUpdate(this.DomainVar, this.ParameterVar, allcomps, extractTaf, SpatialDimension);

            // create derivative operator
            // ==========================

            var JacobianOp = new DifferentialOperator(
                   this.DomainVar,
                   h.JacobianParameterVars,
                   this.CodomainVar,
                   this.QuadOrderFunction);

            if (this.TemporalOperator != null)
                JacobianOp.TemporalOperator = new TemporalOperatorContainer(JacobianOp, this.TemporalOperator);

            foreach (string CodNmn in this.CodomainVar) {
                foreach(var eq in this.EquationComponents[CodNmn]) {

                    if(!(eq is ISupportsJacobianComponent _eq))
                        throw new NotSupportedException(string.Format("Unable to handle component {0}: To obtain a Jacobian operator, all components must implement the {1} interface.", eq.GetType().Name, typeof(ISupportsJacobianComponent).Name));
                    bool eq_suppCoeffUpd = eq is IEquationComponentCoefficient;

                    foreach(var eqj in _eq.GetJacobianComponents(SpatialDimension)) {
                        bool eqj_suppCoeffUpd = eqj is IEquationComponentCoefficient;
                        if(eq_suppCoeffUpd && !eqj_suppCoeffUpd)
                            throw new NotSupportedException("Form '" + eq.GetType().Name + "' supports '" + typeof(IEquationComponentCoefficient).Name + "', but Jacobian Form '" + eqj.GetType().Name + "' does not!");

                        JacobianOp.EquationComponents[CodNmn].Add(eqj);
                    }
                }
            }



            // return
            // =====

            foreach(string domName in this.DomainVar)
                JacobianOp.FreeMeanValue[domName] = this.FreeMeanValue[domName];
            JacobianOp.EdgeQuadraturSchemeProvider = this.EdgeQuadraturSchemeProvider;
            JacobianOp.VolumeQuadraturSchemeProvider = this.VolumeQuadraturSchemeProvider;
            foreach(var kv in this.UserDefinedValues)
                JacobianOp.UserDefinedValues.Add(kv);

            JacobianOp.LinearizationHint = LinearizationHint.AdHoc;

            foreach(DelParameterFactory f in this.ParameterFactories) 
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
            foreach(var cod in this.CodomainVar) {
                foreach(var comp in this.EquationComponents[cod]) {
                    if(comp is IDGdegreeConstraint dgconstr) {
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

        /// <summary>
        /// Used by <see cref="_GetJacobiOperator(int)"/> to encalsulate the temporal operator
        /// of this operator (because of the ownership, the temporal operator cannot be reused).
        /// </summary>
        class TemporalOperatorContainer : ITemporalOperator {

            DifferentialOperator m_newOwner;
            ITemporalOperator m_encapsulatedObj;
            public TemporalOperatorContainer(DifferentialOperator __newOwner, ITemporalOperator __encapsulatedObj) {
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
                if(m_FreeMeanValue == null) {
                    m_FreeMeanValue = new MyDict(this);
                }
                return m_FreeMeanValue;
            }
        }


        /// <summary>
        /// I hate shit like this class - so many dumb lines of code.
        /// </summary>
        class MyDict : IDictionary<string, bool> {

            DifferentialOperator owner;

            public MyDict(DifferentialOperator __owner) {
                owner = __owner;
                InternalRep = new Dictionary<string, bool>();
                foreach(string domName in __owner.DomainVar) {
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
                (InternalRep as ICollection<KeyValuePair<string,bool>>).CopyTo(array, arrayIndex);
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
    }
}
