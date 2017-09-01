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
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;
using BoSSS.Platform;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Stores information about a session.
    /// </summary>
    [Serializable]
    [DataContract]
    public partial class SessionInfo : ISessionInfo {

        /// <summary>
        /// Creates a new instance of SessionInfo.
        /// </summary>
        /// <param name="uid">The unique identifier of the session.</param>
        /// <param name="database">The database where the session is stored.</param>
        public SessionInfo(Guid uid, IDatabaseInfo database) {
            ID = uid;
            Database = database;
            CreationTime = DateTime.Now;

            // initialize tags as empty list
            Tags = new List<string>();
        }

        /// <summary>
        /// Description for the session. Changing this value causes an 
        /// immediate IO operation.
        /// </summary>
        public string Description {
            get {
                return m_Description;
            }
            set {
                m_Description = value == null ? "" : value.Trim();
                this.Save();
            }
        }

        [DataMember]
        private string m_Description;

        /// <summary>
        /// Associates this session with a project.  Changing this value causes
        /// an  immediate IO operation.
        /// </summary>
        public string ProjectName {
            get {
                return m_ProjectName;
            }
            set {
                if (String.IsNullOrWhiteSpace(value)) {
                    throw new Exception("New project name is invalid.");
                }

                m_ProjectName = value.Trim();
                this.Save();
            }
        }

        [DataMember]
        private string m_ProjectName;

        /// <summary>
        /// All the time-steps of this session.
        /// </summary>
        public IList<ITimestepInfo> Timesteps {
            get {
                return Database.Controller.GetTimestepInfos(this).OrderBy(t => t.WriteTime).ToList();
            }
        }

        [DataMember]
        private KeysDict m_KeysAndQueries;

        /// <summary>
        /// see <see cref="ISessionInfo.KeysAndQueries"/>
        /// </summary>
        public IDictionary<string, object> KeysAndQueries {
            get {
                if (m_KeysAndQueries == null) {
                    m_KeysAndQueries = new KeysDict() { m_Owner = this };
                }
                return m_KeysAndQueries;
            }
        }



        ///// <summary>
        ///// A list of all time-steps, initialized lazily
        ///// </summary>
        //[NonSerialized]
        //private IDictionary<Guid, ITimestepInfo> m_Timesteps = new Dictionary<Guid, ITimestepInfo>();

        //[OnDeserialized]
        //private void OnDeserialized(StreamingContext context) {
        //    m_Timesteps = new Dictionary<Guid, ITimestepInfo>();
        //}

        /// <summary>
        /// The session ID this session has been restarted from.
        /// </summary>
        public Guid RestartedFrom {
            get {
                return m_RestartedFrom;
            }
            set {
                m_RestartedFrom = value;
            }
        }

        [DataMember]
        private Guid m_RestartedFrom;

        /// <summary>
        /// A collection of tags for this session.
        /// </summary>
        public IEnumerable<string> Tags {
            get {
                return m_Tags;
            }
            set {
                if (value != null) {
                    m_Tags = value.Where(tag => !String.IsNullOrWhiteSpace(tag))
                        .Distinct().ToArray();
                    this.Save();
                }
            }
        }

        [DataMember]
        private string[] m_Tags;

        /// <summary>
        /// Returns all the grids used in this session
        /// </summary>
        public IEnumerable<IGridInfo> GetGrids() {
            return Database.Controller.GetGridInfos(this);
        }

        #region IDatabaseEntityInfo<ISessionInfo> Members

        /// <summary>
        /// Unique identifier of the session.
        /// </summary>
        public Guid ID {
            get {
                return m_ID;
            }
            private set {
                m_ID = value;
            }
        }

        [DataMember]
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
        /// The name of the session.
        /// </summary>
        public string Name {
            get {
                return m_Name;
            }
            set {
                if (String.IsNullOrWhiteSpace(value)) {
                    throw new Exception("New name of session is invalid.");
                }

                m_Name = value.Trim();
                this.Save();
            }
        }

        [DataMember]
        private string m_Name;

        /// <summary>
        /// The database where this info object is located.
        /// </summary>
        public IDatabaseInfo Database {
            get {
                return m_Database;
            }
            set {
                if (value == null) {
                    m_Database = NullDatabaseInfo.Instance;
                } else {
                    m_Database = value;
                }
            }
        }

        [NonSerialized]
        private IDatabaseInfo m_Database;

        /// <summary>
        /// Copies this IDatabaseObjectInfo object for storage in a different database.
        /// </summary>
        /// <param name="targetDatabase">The target database</param>
        /// <returns>
        /// A copy of this IDatabaseObjectInfo object with all the same
        /// information, except for the database field, which will be the one of the 
        /// target database
        /// </returns>
        public ISessionInfo CopyFor(IDatabaseInfo targetDatabase) {
            SessionInfo copy = new SessionInfo(this.ID, targetDatabase);

            // Copy all the field values
            foreach (var field in this.GetType().GetFields(BindingFlags.Public
                | BindingFlags.NonPublic | BindingFlags.Instance)) {
                field.SetValue(copy, field.GetValue(this));
            }

            copy.m_Database = targetDatabase; // only field with different value
            copy.Save(); // to commit changes

            return copy;
        }

        #endregion

        #region Object Members

        /// <summary>
        /// Creates a string representation of this object.
        /// </summary>
        /// <returns>A string representing the object.</returns>
        public override string ToString() {
            return "{ Guid = " + ID + "; " + CreationTime.ToString() + " Name = " + Name + " }";
        }

        #endregion

        #region IEquatable<ISessionInfo> Members

        /// <summary>
        /// Compares this session to another.
        /// </summary>
        /// <param name="other">The session to compare to.</param>
        /// <returns>
        /// true, if the session GUIDs are the same; false otherwise.
        /// </returns>
        public bool Equals(ISessionInfo other) {
            if (other == null) {
                return false;
            } else {
                return this.ID.Equals(other.ID);
            }
        }

        #endregion

        /// <summary>
        /// Names of compute nodes on which the session is running; Index: MPI
        /// rank index in the MPI_COMM_WORLD communicator;
        /// </summary>
        public IList<string> ComputeNodeNames {
            get {
                return m_ComputeNodeNames;
            }
        }

        /// <summary>
        /// <see cref="ComputeNodeNames"/>
        /// </summary>
        [DataMember]
        private List<string> m_ComputeNodeNames = new List<string>();

        [NonSerialized]
        TextWriter m_TimeStepLog;

        internal void LogTimeStep(Guid g) {
            if (m_TimeStepLog == null) {
                m_TimeStepLog = Database.Controller.DBDriver.FsDriver
                    .GetNewLog("TimestepLog", this.ID);
            }

            m_TimeStepLog.WriteLine(g.ToString());
            m_TimeStepLog.Flush();
        }

        /// <summary>
        /// Saves the current state of this object to the associated
        /// <see cref="Database"/>.
        /// </summary>
        /// <seealso cref="IDatabaseController.SaveSessionInfo"/>
        public void Save() {
            if (Database.Controller.DBDriver.MyRank == 0) {
                if (!this.ID.Equals(Guid.Empty)) {
                    Database.Controller.SaveSessionInfo(this);
                }
            }
        }

        /// <summary>
        /// Dispo.
        /// </summary>
        public void Dispose() {
            if (m_TimeStepLog != null) {
                m_TimeStepLog.Flush();
                m_TimeStepLog.Close();
                m_TimeStepLog.Dispose();
            }
        }
    }
}
