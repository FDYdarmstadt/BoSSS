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
using BoSSS.Platform;
using System.Diagnostics;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// A proxy for <see cref="SessionInfo"/> objects that allows for lazy
    /// loading
    /// </summary>
    public class SessionProxy : ISessionInfo {

        /// <summary>
        /// The real session reflected by this object
        /// </summary>
        private ExpirableLazy<SessionInfo> realSessionInfo;

        /// <summary>
        /// The real session reflected by this object
        /// </summary>
        public ISessionInfo RealSessionInfo {
            get {
                // Allow graceful handling when loading faulty sessions
                ISessionInfo nullSession = NullSessionInfo.Instance;
                return realSessionInfo.Value ?? nullSession;
            }
        }

        /// <summary>
        /// Dissing.
        /// </summary>
        public void Dispose() {
            if(realSessionInfo.IsValueCreated) {
                realSessionInfo.Value.Dispose();
            }
        }


        /// <summary>
        /// Constructs a proxy for the session with id
        /// <paramref name="sessionID"/> within the given
        /// <paramref name="database"/>
        /// </summary>
        /// <param name="sessionID"></param>
        /// <param name="database"></param>
        public SessionProxy(Guid sessionID, IDatabaseInfo database) {
            this.ID = sessionID;
            this.Database = database;
            realSessionInfo = new ExpirableLazy<SessionInfo>(
                delegate () {
                    // Allow graceful handling when loading faulty sessions
                    try {
                        return database.Controller.DBDriver.LoadSession(sessionID, database);
                    } catch(Exception e) {
                        Console.WriteLine(
                            "Loading session {0} failed with message '{1}'", ID, e.Message);
                        return null;
                    }
                },
                s => Utils.GetSessionFileWriteTime(s) == s.WriteTime);
        }

        /// <summary>
        /// See <see cref="SessionInfo.ToString"/>
        /// </summary>
        /// <returns></returns>
        public override string ToString() {
            return RealSessionInfo.ToString();
        }

        #region ISessionInfo Members

        /// <summary>
        /// see <see cref="ISessionInfo.SuccessfulTermination"/>
        /// </summary>
        public bool SuccessfulTermination {
            get {
                return RealSessionInfo.SuccessfulTermination;
            }
        }

        /// <summary>
        /// See <see cref="SessionInfo.Description"/>
        /// </summary>
        public string Description {
            get {
                return RealSessionInfo.Description;
            }
            set {
                RealSessionInfo.Description = value;
            }
        }

        /// <summary>
        /// See <see cref="SessionInfo.ProjectName"/>
        /// </summary>
        public string ProjectName {
            get {
                return RealSessionInfo.ProjectName;
            }
            set {
                RealSessionInfo.ProjectName = value;
            }
        }

        /// <summary>
        /// See <see cref="SessionInfo.Timesteps"/>
        /// </summary>
        public IList<ITimestepInfo> Timesteps {
            get {
                return RealSessionInfo.Timesteps;
            }
        }

        /// <summary>
        /// See <see cref="SessionInfo.RestartedFrom"/>
        /// </summary>
        public Guid RestartedFrom {
            get {
                return RealSessionInfo.RestartedFrom;
            }
        }

        /// <summary>
        /// See <see cref="SessionInfo.MasterGitCommit"/>
        /// </summary>
        public string MasterGitCommit {
            get {
                return RealSessionInfo.MasterGitCommit;
            }
        }

        /// <summary>
        /// See <see cref="SessionInfo.Tags"/>
        /// </summary>
        public IEnumerable<string> Tags {
            get {
                return RealSessionInfo.Tags;
            }
            set {
                RealSessionInfo.Tags = value;
            }
        }

        /// <summary>
        /// See <see cref="SessionInfo.ComputeNodeNames"/>
        /// </summary>
        public IList<string> ComputeNodeNames {
            get {
                return RealSessionInfo.ComputeNodeNames;
            }
        }

        /// <summary>
        /// See <see cref="SessionInfo.GetGrids"/>
        /// </summary>
        public IEnumerable<IGridInfo> GetGrids() {
            return RealSessionInfo.GetGrids();
        }

        #endregion

        #region IDatabaseEntityInfo<ISessionInfo> Members

        /// <summary>
        /// The id of the session
        /// </summary>
        public Guid ID {
            get;
            private set;
        }

        /// <summary>
        /// See <see cref="SessionInfo.CreationTime"/>
        /// </summary>
        public DateTime CreationTime {
            get {
                return RealSessionInfo.CreationTime;
            }
        }

        /// <summary>
        /// See <see cref="SessionInfo.WriteTime"/>. Note that reading this
        /// property does <b>not</b> induce deserialization of the reflected
        /// session.
        /// </summary>
        public DateTime WriteTime {
            get {
                return Utils.GetSessionFileWriteTime(this);
            }
            set {
                RealSessionInfo.WriteTime = value;
            }
        }

        /// <summary>
        /// See <see cref="SessionInfo.Name"/>
        /// </summary>
        public string Name {
            get {
                return RealSessionInfo.Name;
            }
            set {
                RealSessionInfo.Name = value;
            }
        }

        /// <summary>
        /// The associated database
        /// </summary>
        public IDatabaseInfo Database {
            get;
            private set;
        }

        /// <summary>
        /// see <see cref="ISessionInfo.KeysAndQueries"/>
        /// </summary>
        public IDictionary<string, object> KeysAndQueries {
            get {
                return RealSessionInfo.KeysAndQueries;
            }
        }

        /// <summary>
        /// Creates a new session <b>proxy</b> for the given
        /// <paramref name="targetDatabase"/>.
        /// </summary>
        /// <param name="targetDatabase"></param>
        /// <returns></returns>
        public ISessionInfo CopyFor(IDatabaseInfo targetDatabase) {
            return new SessionProxy(ID, targetDatabase);
        }

        #endregion

        #region IEquatable<ISessionInfo> Members

        /// <summary>
        /// See <see cref="SessionInfo.Equals"/>
        /// </summary>
        public bool Equals(ISessionInfo other) {
            return RealSessionInfo.Equals(other);
        }

        #endregion

        /// <summary>
        /// Null object that is used whenever loading a session fails.
        /// </summary>
        public class NullSessionInfo : ISessionInfo {

            /// <summary>
            /// The one and only instance of this class.
            /// </summary>
            public static NullSessionInfo Instance = new NullSessionInfo();

            /// <summary>
            /// Forbid creation
            /// </summary>
            private NullSessionInfo() {
            }

            #region ISessionInfo Members

            /// <summary>
            /// Always empty; writing is a non-operation
            /// </summary>
            public string Description {
                get {
                    return "";
                }
                set {
                    //nop
                }
            }

            /// <summary>
            /// Always empty; writing is a non-operation
            /// </summary>
            public string ProjectName {
                get {
                    return "";
                }
                set {
                    //nop
                }
            }

            /// <summary>
            /// An empty list
            /// </summary>
            public IList<ITimestepInfo> Timesteps {
                get {
                    return new List<ITimestepInfo>();
                }
            }

            /// <summary>
            /// An empty guid
            /// </summary>
            public Guid RestartedFrom {
                get {
                    return Guid.Empty;
                }
            }

            /// <summary>
            /// An empty string
            /// </summary>
            public string MasterGitCommit {
                get {
                    return "";
                }
            }

            /// <summary>
            /// 
            /// </summary>
            public IEnumerable<string> Tags {
                get {
                    yield break;
                }
                set {
                    //nop
                }
            }

            /// <summary>
            /// An empty list
            /// </summary>
            public IList<string> ComputeNodeNames {
                get {
                    return new List<string>();
                }
            }

            /// <summary>
            /// An empty list
            /// </summary>
            /// <returns></returns>
            public IEnumerable<IGridInfo> GetGrids() {
                yield break;
            }

            #endregion

            #region IDatabaseEntityInfo<ISessionInfo> Members

            /// <summary>
            /// An empty guid
            /// </summary>
            public Guid ID {
                get {
                    return Guid.Empty;
                }
            }

            /// <summary>
            /// See <see cref="DateTime.MinValue"/>
            /// </summary>
            public DateTime CreationTime {
                get {
                    return DateTime.MinValue;
                }
            }

            /// <summary>
            /// See <see cref="DateTime.MinValue"/>
            /// </summary>
            public DateTime WriteTime {
                get {
                    return DateTime.MinValue;
                }
                set {
                    throw new NotImplementedException();
                }
            }

            /// <summary>
            /// An empty string; writing is a non-operation
            /// </summary>
            public string Name {
                get {
                    return "";
                }
                set {
                    //nop
                }
            }

            /// <summary>
            /// See <see cref="NullDatabaseInfo.Instance"/>
            /// </summary>
            public IDatabaseInfo Database {
                get {
                    return NullDatabaseInfo.Instance;
                }
            }

            /// <summary>
            /// Nix
            /// </summary>
            public IDictionary<string, object> KeysAndQueries {
                get {
                    return new SortedDictionary<string, object>();
                }
            }

            /// <summary>
            /// always false
            /// </summary>
            public bool SuccessfulTermination {
                get {
                    return false;
                }
            }

            /// <summary>
            /// Always throws <see cref="NotImplementedException"/>
            /// </summary>
            /// <param name="targetDatabase"></param>
            /// <returns></returns>
            public ISessionInfo CopyFor(IDatabaseInfo targetDatabase) {
                throw new NotImplementedException();
            }

            #endregion

            #region IEquatable<ISessionInfo> Members

            /// <summary>
            /// See <see cref="object.ReferenceEquals"/>
            /// </summary>
            /// <param name="other"></param>
            /// <returns></returns>
            public bool Equals(ISessionInfo other) {
                return object.ReferenceEquals(this, other);
            }

            #endregion

            /// <summary>
            /// Displays an error message
            /// </summary>
            /// <returns></returns>
            public override string ToString() {
                return "{ Error loading session }";
            }

            /// <summary>
            /// Does nothing.
            /// </summary>
            public void Dispose() {
            }
        }
    }
}
