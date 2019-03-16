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
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;
using BoSSS.Foundation.Grid;
using ilPSP;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Information about a single time-step
    /// </summary>
    [Serializable]
    [DataContract]
    public class TimestepInfo : ITimestepInfo {

        #region Constructors

        /// <summary>
        /// empty constructor for serialization
        /// </summary>
        protected TimestepInfo() {
        }

        /// <summary>
        /// Constructs information about a time-step with given
        /// <paramref name="fields"/>.
        /// </summary>
        /// <param name="physTime">
        /// The physical time represented by this object
        /// </param>
        /// <param name="session">
        /// The session this time-step belongs to
        /// </param>
        /// <param name="TimestepNo">
        /// The number of the represented time-step
        /// </param>
        /// <param name="fields">
        /// The fields associated with this time-step.
        /// </param>
        /// <param name="storageID">
        /// <see cref="Guid"/> of the storage vector containing the serialized
        /// information stored in this object.
        /// </param>
        public TimestepInfo(
            double physTime, ISessionInfo session, TimestepNumber TimestepNo, IEnumerable<DGField> fields) {

            // check & set grid
            this.m_GridGuid = fields.Count() > 0 ? fields.First().Basis.GridDat.GridID : Guid.Empty;
            if(fields.Where(f => !f.Basis.GridDat.GridID.Equals(this.m_GridGuid)).Count() > 0) {
                throw new ArgumentException("all fields must be associated to the same grid.");
            }

            // check validity of identifications: 
            {
                List<DGField> _Fields = new List<DGField>();
                FlattenHierarchy(_Fields, fields);
                HashSet<string> allIdentifications = new HashSet<string>(new ilPSP.FuncEqualityComparer<string>((a, b) => a.Equals(b), a => a.GetHashCode()));
                foreach (var f in _Fields) {
                    if (f.Identification == null || f.Identification.Length <= 0) {
                        throw new ArgumentException("Creating timesteps from fields without an identification is not supported.");
                    }

                    if (allIdentifications.Contains(f.Identification)) {
                        throw new ArgumentException("Within a Timestep, the Identification of each field must be unique: found at least two fields named '" + f.Identification + "'");
                    }

                    allIdentifications.Add(f.Identification);
                }
            }

            // set members:
            ID = Guid.Empty;
            this.m_TimestepNumber = TimestepNo;
            this.PhysicalTime = physTime;
            this.Session = session;
            this.m_FieldInitializers = fields.Select(f => f.Initializer).ToArray();
            CreationTime = DateTime.Now;
            this.m_StorageID = Guid.Empty;
        }

        /// <summary>
        /// Creates a new instance of <see cref="TimestepInfo"/>
        /// </summary>
        /// <param name="uid">The unique identifier of this time-step.</param>
        /// <param name="database">
        /// The database this time-step is associated with.
        /// </param>
        public TimestepInfo(Guid uid, IDatabaseInfo database) {
            ID = uid;
            m_Database = database;
        }

        #endregion

        /// <summary>
        /// The 'flattening' of the hierarchy is required because XDG fields
        /// depend on Level Sets.
        /// </summary>
        internal static void FlattenHierarchy(List<DGField> output, IEnumerable<DGField> input) {
            foreach (DGField f in input) {
                // recursion
                FlattenHierarchy(output, f.ReportDependentFields());
                if (!output.Contains(f, (a, b) => object.ReferenceEquals(a, b))) {
                    output.Add(f);
                }
            }
        }

        /// <summary>
        /// Initializer of the fields associated with this time-step
        /// </summary>
        public IEnumerable<DGField.FieldInitializer> FieldInitializers {
            get {
                return this.m_FieldInitializers;
            }
        }

        [DataMember]
        private DGField.FieldInitializer[] m_FieldInitializers;

        /// <summary>
        /// Guid of the vector which contains all the data of the time-step
        /// - empty after construction
        /// - set through <see cref="DatabaseDriver.SaveTimestep"/>
        /// </summary>
        public Guid StorageID {
            get {
                return m_StorageID;
            }
            internal set {
                m_StorageID = value;
            }
        }

        [DataMember]
        private Guid m_StorageID;

        /// <summary>
        /// See <see cref="ITimestepInfo"/>
        /// </summary>
        public TimestepNumber TimeStepNumber {
            get {
                return this.m_TimestepNumber;
            }
        }

        [DataMember]
        private TimestepNumber m_TimestepNumber;

        /// <summary>
        /// The grid associated with this time-step
        /// </summary>
        public IGridInfo Grid {
            get {
                if (Database != null) {
                    return Database.Controller.GetGridInfo(this.m_GridGuid);
                } else {
                    throw new NullReferenceException("Database for grid " +
                        this.ID + " not set.");
                }
            }
        }

        /// <summary>
        /// ID of the grid associated with this time-step
        /// </summary>
        public Guid GridID {
            get {
                return m_GridGuid;
            }
        }

        [DataMember]
        private Guid m_GridGuid;

        [NonSerialized]
        private ISessionInfo session;

        /// <summary>
        /// The associated session
        /// </summary>
        public ISessionInfo Session {
            get {
                return session;
            }
            set {
                session = value;
            }
        }

        /// <summary>
        /// Physical time of the time-step
        /// </summary>
        public double PhysicalTime {
            get {
                return m_PhysicalTime;
            }
            private set {
                m_PhysicalTime = value;
            }
        }

        [DataMember]
        private double m_PhysicalTime;

        /// <summary>
        /// Contains information on the fields of this time-step.
        /// </summary>
        public IEnumerable<DGField> Fields {
            get {
                return m_Fields.Value;
            }
        }

        [NonSerialized]
        private Lazy<IEnumerable<DGField>> m_Fields;

        /// <summary>
        /// Initialization method that is used both by the constructors
        /// and by the serialization to initialize non-persistent fields.
        /// </summary>
        /// <param name="context"></param>
        [OnDeserialized]
        private void Initialize(StreamingContext context) {
            // init fields as lazy list
            m_Fields = new Lazy<IEnumerable<DGField>>(GetFields);
        }

        /// <summary>
        /// Populates the lazy list of fields.
        /// </summary>
        /// <returns></returns>
        private IEnumerable<DGField> GetFields() {
            IInitializationContext ic = this.Database.Controller.GetInitializationContext(this);

            List<DGField> fields = new List<DGField>();
            foreach (DGField.FieldInitializer fi in this.m_FieldInitializers) {
                DGField fld = fi.Initialize(ic);
                fields.Add(fld);
            }
            Database.Controller.DBDriver.LoadFieldData(this, ic.GridData, fields);
            return fields;
        }

        #region IDatabaseEntityInfo<ITimestepInfo> Members

        /// <summary>
        /// Unique identifier of the TimestepInfo object.
        /// - empty after construction
        /// - set through <see cref="DatabaseDriver.SaveTimestep"/>
        /// </summary>
        public Guid ID {
            get {
                return m_ID;
            }
            internal set {
                m_ID = value;
            }
        }

        [NonSerialized]
        private Guid m_ID;

        /// <summary>
        /// The time when the represented entity has been created.
        /// </summary>
        public DateTime CreationTime {
            get {
                return m_CreationTime;
            }
            private set {
                m_CreationTime = value;
            }
        }

        [DataMember]
        private DateTime m_CreationTime;

        /// <summary>
        /// The time when this object has been written to disc.
        /// </summary>
        public DateTime WriteTime {
            get {
                return m_WriteTime;
            }
            set {
                m_WriteTime = value;
            }
        }

        [NonSerialized]
        DateTime m_WriteTime;

        /// <summary>
        /// The name of the entity.
        /// </summary>
        public string Name {
            get {
                return m_Name;
            }
            set {
                if (value == null) {
                    throw new ArgumentNullException();
                }

                m_Name = value.Trim(); // allow time-step name to be empty
            }
        }

        [DataMember]
        private string m_Name;

        /// <summary>
        /// The database where this info object is located
        /// </summary>
        public IDatabaseInfo Database {
            get {
                return m_Database;
            }
            set {
                if (value != null) {
                    m_Database = value;
                } else {
                    throw new ArgumentNullException();
                }
            }
        }

        [NonSerialized]
        private IDatabaseInfo m_Database;

        /// <summary>
        /// Copies this ITimestepInfo object for storage in a different database.
        /// </summary>
        /// <param name="targetDatabase">The target database</param>
        /// <returns>
        /// A copy of this ITimestepInfo object with all the same
        /// information, except for the database field, which will be the one of the 
        /// target database
        /// </returns>
        public ITimestepInfo CopyFor(IDatabaseInfo targetDatabase) {
            TimestepInfo copy = new TimestepInfo(ID, Database);

            // Copy all the field values
            foreach (var field in this.GetType().GetFields(BindingFlags.Public
                | BindingFlags.NonPublic | BindingFlags.Instance)) {
                field.SetValue(copy, field.GetValue(this));
            }

            copy.m_Database = targetDatabase;
            copy.m_TimestepNumber = this.TimeStepNumber;

            return copy;
        }

        #endregion

        /// <summary>
        /// Reports time-step number and physical time
        /// </summary>
        /// <returns></returns>
        public override string ToString() {
            return String.Format(
                " {{ Time-step: {0}; Physical time: {1}s; Fields: {2}; Name: {3} }}",
                TimeStepNumber.ToString(),
                PhysicalTime,
                FieldInitializers.Select(f => f.Identification).Aggregate((s, t) => s + ", " + t),
                Name);
        }
    }
}
