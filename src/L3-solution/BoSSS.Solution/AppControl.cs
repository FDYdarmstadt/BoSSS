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
using BoSSS.Foundation.IO;
using BoSSS.Solution.Queries;
using ilPSP;
using System.Linq;
using System.Reflection;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;

namespace BoSSS.Solution.Control {

    /// <summary>
    /// Version 2 of the Application control
    /// </summary>
    [Serializable]
    public class AppControl {


        public abstract Application<T> GetSolver<T>() where T : AppControl;

        /// <summary>
        /// The generating code as text representation
        /// </summary>
        public string ControlFileText {
            get;
            internal set;
        }

        /// <summary>
        /// ctor.
        /// </summary>
        public AppControl() {
            this.FieldOptions = new Dictionary<string, FieldOpts>();
            this.BoundaryValues = new Dictionary<string, BoundaryValueCollection>(new InvariantCultureIgnoreCase_Comparer());
            this.Tags = new List<string>();
            this.m_InitialValues_Evaluators = new Dictionary<string, Func<double[], double>>();
            this.m_InitialValues = new Dictionary<string, IBoundaryAndInitialData>();
            this.NoOfMultigridLevels = 0;
        }

        [Serializable]
        class InvariantCultureIgnoreCase_Comparer : IEqualityComparer<string> {
            public bool Equals(string a, string b) {
                return a.Equals(b, StringComparison.InvariantCultureIgnoreCase);
            }

            public int GetHashCode(string obj) {
                return obj.GetHashCode();
            }
        }


        /// <summary>
        /// DG polynomial degree for fields;
        ///  - key: string, the field identification string
        ///  - value: options for the DG field
        /// </summary>
        public IDictionary<string, FieldOpts> FieldOptions {
            get;
            private set;
        }

        /// <summary>
        /// Adds an entry to <see cref="FieldOptions"/>.
        /// </summary>
        /// <param name="Degree">
        /// Polynomial degree of ther DG field; if negative, the application tries to determine a degree by itself.
        /// </param>
        public void AddFieldOption(string DGFieldName, int Degree = -1, FieldOpts.SaveToDBOpt SaveOpt = FieldOpts.SaveToDBOpt.TRUE) {
            FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = Degree,
                SaveToDB = SaveOpt
            });
        }

        /// <summary>
        /// Checks that all control file options satisfy the conditions
        /// specified in the corresponding
        /// <see cref="ControlOptionRequirementAttribute"/>s. Override this
        /// method to implement more complicated checks that e.g. depend on the
        /// values of other attributes. For example, one could verify that an
        /// LBB-condition is fulfilled.
        /// </summary>
        public virtual void Verify() {
            Type currentType = this.GetType();
            List<Type> types = new List<Type>() { currentType };
            while (currentType.BaseType != null) {
                currentType = currentType.BaseType;
                types.Add(currentType);
            }

            foreach (Type type in types) {
                foreach (FieldInfo field in type.GetFields()) {
                    object value = field.GetValue(this);
                    foreach (var attribute in field.GetCustomAttributes<ControlOptionRequirementAttribute>()) {
                        attribute.Verify(field.Name, value);
                    }
                }
            }
        }

        /// <summary>
        /// A collection of all boundary values for a specific edge tag,
        /// see <see cref="Evaluators"/>
        /// </summary>
        [Serializable]
        sealed public class BoundaryValueCollection {

            /// <summary>
            /// optional string to encode the type of the boundary
            /// condition (Neumann or Dirichlet b.c., ...)
            /// </summary>
            public string type;

            /// <summary>
            /// ctor.
            /// </summary>
            public BoundaryValueCollection() {
                this.m_BoundaryValues_Evaluators = new Dictionary<string, Func<double[], double, double>>();
                this.m_BoundaryValues = new Dictionary<string, IBoundaryAndInitialData>();
            }

            /// <summary>
            /// Unfortunately, we cannot serialize delegates.
            /// </summary>
            [NonSerialized]
            Dictionary<string, Func<double[], double, double>> m_BoundaryValues_Evaluators;

            /// <summary>
            /// Serializeable representation of boundary data.
            /// </summary>
            Dictionary<string, IBoundaryAndInitialData> m_BoundaryValues;

            /// <summary>
            /// Steady-state or instationary boundary values.
            /// - key: a name for the boundary value, e.g. 'VelocityX'<br/>
            /// - value: some function that maps a space coordinate to some function value, i.e.  \f$ (\vec{x},t) \mapsto f(\vec{x},t)\f$ 
            /// </summary>
            /// <remarks>
            /// Adding delegates directly to this dictionary is possible for backward compatibility reasons,
            /// although this limits some functionality - e.g. the control object is usually not serializeable anymore.
            /// </remarks>
            public IDictionary<string, Func<double[], double, double>> Evaluators {
                get {
                    if (m_BoundaryValues_Evaluators == null)
                        m_BoundaryValues_Evaluators = new Dictionary<string, Func<double[], double, double>>();

                    foreach (string name in m_BoundaryValues.Keys) {
                        if (!m_BoundaryValues_Evaluators.ContainsKey(name)) {
                            m_BoundaryValues_Evaluators.Add(name, m_BoundaryValues[name].Evaluate);
                        }
                    }
                    return m_BoundaryValues_Evaluators;
                }
            }

            /// <summary>
            /// Steady-state or instationary boundary values.
            /// - key: a name for the boundary value, e.g. 'VelocityX'<br/>
            /// - value: some function that maps a space coordinate to some function value, i.e.  \f$ (\vec{x},t) \mapsto f(\vec{x},t)\f$ 
            /// </summary>
            public IDictionary<string, IBoundaryAndInitialData> Values {
                get {
                    return m_BoundaryValues;
                }
            }
        }


        /// <summary>
        /// A mapping from edge tag names to a collection of boundary values.<br/>
        /// key: edge tag names <see cref="GridCommons.EdgeTagNames"/><br/>
        /// value: boundary values for various fields.
        /// </summary>
        public IDictionary<string, BoundaryValueCollection> BoundaryValues {
            get;
            private set;
        }

        /// <summary>
        /// Creates an empty boundary condition
        /// </summary>
        /// <param name="EdgeTagName">Name of the boundary condition</param>
        public void AddBoundaryCondition(string EdgeTagName) {
            if (!this.BoundaryValues.ContainsKey(EdgeTagName))
                this.BoundaryValues.Add(EdgeTagName, new BoundaryValueCollection());

        }

        /// <summary>
        /// Adds a time-dependent boundary condition.
        /// </summary>
        /// <param name="EdgeTagName">Name of the boundary condition</param>
        /// <param name="fieldname">Name of the field for which the boundary condition is valid</param>
        /// <param name="value">Function of the boundary condition</param>
        public void AddBoundaryCondition(string EdgeTagName, string fieldname, Func<double[], double, double> value) {
            if (!this.BoundaryValues.ContainsKey(EdgeTagName))
                this.BoundaryValues.Add(EdgeTagName, new BoundaryValueCollection());

            if (this.BoundaryValues[EdgeTagName].Evaluators.ContainsKey(fieldname)
                //||  this.BoundaryValues[EdgeTagName].TimedepBoundaryValues.ContainsKey(fieldname) 
                )
                throw new ArgumentException(string.Format("Boundary condition for field '{0}' and edge tag name '{1}' already specified.", EdgeTagName, fieldname));

            this.BoundaryValues[EdgeTagName].Evaluators.Add(fieldname, value);
        }

        /// <summary>
        /// Adds a boundary condition that does not depend on time.
        /// </summary>
        /// <param name="EdgeTagName">Name of the boundary condition</param>
        /// <param name="fieldname">Name of the field for which the boundary condition is valid</param>
        /// <param name="value">Function of the boundary condition</param>
        public void AddBoundaryCondition(string EdgeTagName, string fieldname, Func<double[], double> value) {
            AddBoundaryCondition(EdgeTagName, fieldname, (X, t) => value(X));
        }

        /// <summary>
        /// Adds a time-dependent boundary condition.
        /// </summary>
        /// <param name="EdgeTagName">Name of the boundary condition</param>
        /// <param name="fieldname">Name of the field for which the boundary condition is valid</param>
        /// <param name="FormulaText">
        /// Text representation of a delegate, see <see cref="Formula"/>.
        /// </param>
        /// <param name="TimeDependent">
        /// Whether the formula in <paramref name="FormulaText"/> is time-dependent or not, see <see cref="Formula"/>.
        /// </param>
        public void AddBoundaryCondition(string EdgeTagName, string fieldname, string FormulaText, bool TimeDependent) {
            if (!this.BoundaryValues.ContainsKey(EdgeTagName))
                this.BoundaryValues.Add(EdgeTagName, new BoundaryValueCollection());

            if (this.BoundaryValues[EdgeTagName].Evaluators.ContainsKey(fieldname))
                throw new ArgumentException(string.Format("Boundary condition for field '{0}' and edge tag name '{1}' already specified.", EdgeTagName, fieldname));

            this.BoundaryValues[EdgeTagName].Values.Add(fieldname, new Formula(FormulaText, TimeDependent));
        }

        [NonSerialized]
        Dictionary<string, Func<double[], double>> m_InitialValues_Evaluators;

        /// <summary>
        /// Initial values which are set at the start of the simulation.
        /// - key: identification of the DG field, see <see cref="BoSSS.Foundation.DGField.Identification"/>.
        /// - value: data or mathematical expression which is projected at application startup.
        /// </summary>
        /// <remarks>
        /// Adding delegates directly to this dictionary is possible for backward compatibility reasons,
        /// although this limits some functionality - e.g. the control object is usually not serializeable anymore.
        /// </remarks>
        public IDictionary<string, Func<double[], double>> InitialValues_Evaluators {
            get {
                if (m_InitialValues_Evaluators == null)
                    m_InitialValues_Evaluators = new Dictionary<string, Func<double[], double>>();

                foreach (string name in InitialValues.Keys) {
                    if (!m_InitialValues_Evaluators.ContainsKey(name)) {
                        m_InitialValues_Evaluators.Add(name, X => InitialValues[name].Evaluate(X, 0.0));
                    }
                }
                return m_InitialValues_Evaluators;
            }
        }

        Dictionary<string, IBoundaryAndInitialData> m_InitialValues;

        /// <summary>
        /// Initial values which are set at the start of the simulation.
        /// - key: identification of the DG field, see <see cref="BoSSS.Foundation.DGField.Identification"/>.
        /// - value: data or mathematical expression which is projected at application startup.
        /// </summary>
        public IDictionary<string, IBoundaryAndInitialData> InitialValues {
            get {
                return m_InitialValues;
            }
        }


        ///// <summary>
        ///// Adds a time-dependent initial value (e.g. some external force field which may change over time).
        ///// </summary>
        ///// <param name="fieldname">Name of the DG field which should hold the initial value.</param>
        ///// <param name="value">Function of the boundary condition</param>
        //public void AddInitialValue(string fieldname, Func<double[], double> value) {
        //    this.InitialValues_Evaluators.Add(fieldname, value);
        //}

        /// <summary>
        /// Contains a list of queries that shall be evaluated at the end of
        /// the simulation
        /// </summary>
        [NonSerialized]
        public readonly IDictionary<string, Query> Queries = new Dictionary<string, Query>();

        /// <summary>
        /// Saves restart Information: GUID of the restarted session and the time-step
        /// If empty, no restart is done.
        /// </summary>
        public Tuple<Guid, TimestepNumber> RestartInfo;

        /// <summary>
        /// The GUID of the grid to load. 
        /// The use of this member is exclusive with <see cref="GridFunc"/>, i.e. 
        /// this member may only be set unequal to null if <see cref="GridFunc"/> is null.
        /// </summary>
        public Guid GridGuid;

        /// <summary>
        /// Returns a grid object that should be used for the simulation.
        /// The use of this member is exclusive with <see cref="GridGuid"/>, i.e. 
        /// this member may only be set unequal to null if <see cref="GridGuid"/> is null.
        /// </summary>
        [NonSerialized]
        public Func<GridCommons> GridFunc;

        /// <summary>
        /// Algorithm for grid partitioning.
        /// </summary>
        public GridPartType GridPartType = GridPartType.METIS;

        /// <summary>
        /// grid partitioning options
        /// </summary>
        public string GridPartOptions;

        /// <summary>
        /// tags/keywords to describe a solver run, to aid sorting/search for specific 
        /// solver runs in the database.
        /// </summary>
        public ICollection<string> Tags {
            get;
            private set;
        }

        /// <summary>
        /// Optional session name (for better identification).
        /// </summary>
        public string SessionName;

        /// <summary>
        /// Mandatory project name.
        /// </summary>
        public string ProjectName;

        /// <summary>
        /// optional project description.
        /// </summary>
        public string ProjectDescription;

        /// <summary>
        /// number of time-steps that the solver should perform;
        /// A negative value indicates that this value is not set within the control file;
        /// </summary>
        public int NoOfTimesteps = -1;

        /// <summary>
        /// physical time at which the solver terminates;
        /// </summary>
        public double Endtime = double.MaxValue;

        /// <summary>
        /// interval in which "restart files" should be written to the database
        /// </summary>
        public int saveperiod = 1;

        /// <summary>
        /// lower threshold for the time-step
        /// </summary>
        /// <remarks>
        /// A negative value indicates that this is not initialized;
        /// </remarks>
        public double dtMin = -1;

        /// <summary>
        /// upper threshold for the time-step;
        /// A negative value indicates that this is not initialized;
        /// </summary>
        /// <remarks>
        /// A negative value indicates that this is not initialized;
        /// </remarks>
        public double dtMax = -1;

        /// <summary>
        /// Sets/Gets a fixed time-step size.
        /// </summary>
        public double dtFixed {
            get {
                if (dtMin != dtMax) {
                    return double.NaN;
                }
                return dtMin;
            }
            set {
                dtMin = value;
                dtMax = value;
            }
        }

        /// <summary>
        /// Checks if a fixed time-step size has been set and returns this value.
        /// </summary>
        public double GetFixedTimestep() {
            if (dtMin != dtMax) {
                throw new ApplicationException("Fixed time-step required; minimum and maximum time-step size must be set to same value.");
            }
            return dtMin;
        }


        /// <summary>
        /// See <see cref="CompMode"/>.
        /// </summary>
        public enum _CompMode {

            /// <summary>
            /// Instationary/Transient simulation.
            /// </summary>
            Transient,

            /// <summary>
            /// Steady-State calculation.
            /// </summary>
            Steady
        }


        /// <summary>
        /// For solvers which support both, stationary as well as transient simulations, the corresponding switch.
        /// </summary>
        public _CompMode CompMode;


        /// <summary>
        /// true if information should be written to the database, false
        /// if "passive io" (only reading grids, <see cref="BoSSS.Foundation.IO.IFileSystemDriver"/>)
        /// should be used;
        /// </summary>
        public bool savetodb = true;

        /// <summary>
        /// minimum tracer level priority that is written to the trace file
        /// </summary>
        public string TracingNamespaces = null;

        /// <summary>
        /// File system path to database.
        /// </summary>
        public string DbPath = null;


        /// <summary>
        /// location where the result files of a parameter study should be saved
        /// </summary>
        public string logFileDirectory = ".";

        /// <summary>
        /// If a parameter study is run, this member helps to distinct different cases (aka. sessions or runs) 
        /// of the study, which should be described by an enumeration of
        ///  parameter - name/value - pairs.
        /// </summary>
        public IEnumerable<Tuple<string, object>> Paramstudy_CaseIdentification = null;

        /// <summary>
        /// Continue parameter study if one case (aka. sessions or run) throws an exception.
        /// </summary>
        public bool Paramstudy_ContinueOnError = true;

        /// <summary>
        /// Number of aggregation multi-grid levels, <see cref="Application{T}.MultigridLevels"/>.
        /// </summary>
        public int NoOfMultigridLevels {
            get;
            set;
        }

        /// <summary>
        /// If true, a redistribution will be attempted BEFORE the first
        /// time-step starts
        /// </summary>
        public bool DynamicLoadBalancing_RedistributeAtStartup = false;

        /// <summary>
        /// A method that creates a new estimator for the runtime cost of individual cells
        /// </summary>
        public List<Func<IApplication<AppControl>, int, ICellCostEstimator>> DynamicLoadBalancing_CellCostEstimatorFactories =
            new List<Func<IApplication<AppControl>, int, ICellCostEstimator>>();

        /// <summary>
        /// Number of time-steps, after which dynamic load balancing is performed; if negative, dynamic load balancing is turned off.
        /// </summary>
        public int DynamicLoadBalancing_Period = -1;

        /// <summary>
        /// Relative value, which is compared against the relative computation imbalance
        /// \f[ 
        /// \frac{\text{maximum runtime difference over all MPI processors}}{\text{maximum runtime over all MPI processors}} .
        /// \f]
        /// Dynamic load balancing is suppressed if the relative computation imbalance is below this value.
        /// </summary>
        public double DynamicLoadBalancing_ImbalanceThreshold = 0.12;


        /// <summary>
        /// switch for activating adaptive mesh refinement
        /// </summary>
        public bool AdaptiveMeshRefinement = false;


        /// <summary>
        /// Actual type of cut cell quadrature to use; If no XDG, is used, resp. no cut cells are present,
        /// this setting has no effect.
        /// </summary>
        public XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
    }
}
