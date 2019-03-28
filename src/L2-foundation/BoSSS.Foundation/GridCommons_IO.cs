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
using System.Reflection;
using System.Runtime.Serialization;
using BoSSS.Foundation.IO;
using Newtonsoft.Json;

namespace BoSSS.Foundation.Grid.Classic {

    public partial class GridCommons : IGridInfo, ICloneable, IEquatable<IGridInfo> {

        /// <summary>
        /// creation time of the grid: implements <see cref="IDatabaseEntityInfo{T}.CreationTime"/>
        /// </summary>
        public DateTime CreationTime {
            get {
                return m_CreationTime;
            }
        }

        [DataMember]
        DateTime m_CreationTime;

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

        #region IO.IGridInfo members

        /// <summary>
        /// see <see cref="ID"/>;
        /// </summary>
        [DataMember]
        Guid m_GridGuid;

        /// <summary>
        /// Guid/Identification of this grid object in the database <see cref="Database"/>
        /// </summary>
        public Guid ID {
            get {
                return this.m_GridGuid;
            }
        }

        /// <summary>
        /// grid name: implementation of <see cref="IDatabaseEntityInfo{T}.Name"/>
        /// </summary>
        public string Name {
            get {
                return m_Name;
            }
            set {
                if (String.IsNullOrWhiteSpace(value) == false) {
                    m_Name = value.Trim();
                    if (Database != null) {
                        Database.Controller.SaveGridInfo(this);
                    }
                } else {
                    throw new Exception("New name of grid is invalid.");
                }
            }
        }

        [DataMember]
        private string m_Name;

        /// <summary>
        /// number of cells in the grid: implementation of <see cref="IGridInfo.NumberOfCells"/>
        /// </summary>
        public int NumberOfCells {
            get {
                return (int)NumberOfCells_l;
            }
        }

        [NonSerialized]
        internal IO.IDatabaseInfo m_Database = null;

        /// <summary>
        /// Database which de-serialized this grid; implementation of <see cref="IDatabaseEntityInfo{t}.Database"/>
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

        /// <summary>
        /// Copies this info object for usage in another database.
        /// </summary>
        /// <param name="targetDatabase">The target database.</param>
        /// <returns>
        /// A copy of the original info object for usage in the target database.
        /// </returns>
        public IGridInfo CopyFor(IDatabaseInfo targetDatabase) {
            GridCommons copy = new GridCommons(m_RefElements, m_EdgeRefElements);

            Type type = this.GetType();
            // Cycle to get all the fields from the parent classes as well.
            while (type != null) {
                // Copy all the field values
                foreach (var field in type.GetFields(BindingFlags.Public
                    | BindingFlags.NonPublic | BindingFlags.Instance)) {
                    field.SetValue(copy, field.GetValue(this));
                }
                type = type.BaseType;
            }
            copy.m_Database = targetDatabase;

            return copy;
        }

        [NonSerialized]
        [JsonIgnore]
        GridCommonsDatabaseMethods dataBaseMethods;

        public IGridSerializationHandler GridSerializationHandler {
            get {
                if (dataBaseMethods == null)
                    dataBaseMethods = new GridCommonsDatabaseMethods(this);
                return dataBaseMethods;
            }
        }

        #endregion

        #region Object members
        /// <summary>
        /// Creates a human-readable string representation of the grid info object.
        /// </summary>
        /// <returns>A string summary of the grid.</returns>
        public override string ToString() {
            return "{ Guid = " + ID + "; Name = " + Name + "; Cell Count = "
                + NumberOfCells + "; Dim = " + SpatialDimension + " }";
        }

        #endregion
    }
}
