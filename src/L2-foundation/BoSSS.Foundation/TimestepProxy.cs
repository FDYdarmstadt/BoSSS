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

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// A proxy for <see cref="TimestepInfo"/> objects that allows for lazy
    /// loading
    /// </summary>
    public class TimestepProxy : ITimestepInfo {

        /// <summary>
        /// The real time-step reflected by this object
        /// </summary>
        private ExpirableLazy<TimestepInfo> realTimestepInfo;

        /// <summary>
        /// Constructs a proxy for the time-step with the given
        /// <paramref name="timestepGuid"/> within the given
        /// <paramref name="session"/>.
        /// </summary>
        /// <param name="timestepGuid"></param>
        /// <param name="session"></param>
        public TimestepProxy(Guid timestepGuid, ISessionInfo session) {
            this.ID = timestepGuid;
            this.Session = session;
            realTimestepInfo = new ExpirableLazy<TimestepInfo>(
                () => session.Database.Controller.DBDriver.LoadTimestepInfo(timestepGuid, session, session.Database),
                t => Utils.GetTimestepFileWriteTime(t) == t.WriteTime);
        }

        /// <summary>
        /// See <see cref="TimestepInfo.ToString"/>
        /// </summary>
        /// <returns></returns>
        public override string ToString() {
            return realTimestepInfo.Value.ToString();
        }

        #region ITimestepInfo Members

        /// <summary>
        /// See <see cref="TimestepInfo.TimeStepNumber"/>
        /// </summary>
        public TimestepNumber TimeStepNumber {
            get {
                return realTimestepInfo.Value.TimeStepNumber;
            }
        }

        /// <summary>
        /// See <see cref="TimestepInfo.Grid"/>
        /// </summary>
        public IGridInfo Grid {
            get {
                return realTimestepInfo.Value.Grid;
            }
        }

        /// <summary>
        /// The associated session
        /// </summary>
        public ISessionInfo Session {
            get;
            private set;
        }

        /// <summary>
        /// See <see cref="TimestepInfo.GridID"/>
        /// </summary>
        public Guid GridID {
            get {
                return realTimestepInfo.Value.GridID;
            }
        }

        /// <summary>
        /// See <see cref="TimestepInfo.PhysicalTime"/>
        /// </summary>
        public double PhysicalTime {
            get {
                return realTimestepInfo.Value.PhysicalTime;
            }
        }

        /// <summary>
        /// See <see cref="TimestepInfo.StorageID"/>
        /// </summary>
        public Guid StorageID {
            get {
                return realTimestepInfo.Value.StorageID;
            }
        }

        /// <summary>
        /// See <see cref="TimestepInfo.FieldInitializers"/>
        /// </summary>
        public IEnumerable<DGField.FieldInitializer> FieldInitializers {
            get {
                return realTimestepInfo.Value.FieldInitializers;
            }
        }

        /// <summary>
        /// See <see cref="TimestepInfo.Fields"/>
        /// </summary>
        public IEnumerable<DGField> Fields {
            get {
                return realTimestepInfo.Value.Fields;
            }
        }

        #endregion

        #region IDatabaseEntityInfo<ITimestepInfo> Members

        /// <summary>
        /// The unique ID of the time-step
        /// </summary>
        public Guid ID {
            get;
            private set;
        }

        /// <summary>
        /// See <see cref="TimestepInfo.CreationTime"/>
        /// </summary>
        public DateTime CreationTime {
            get {
                return realTimestepInfo.Value.CreationTime;
            }
        }

        /// <summary>
        /// See <see cref="TimestepInfo.WriteTime"/>
        /// </summary>
        /// <remarks>
        /// Note that retrieving this value does not cause deserialization of
        /// <see cref="realTimestepInfo"/>, while setting this value obviously
        /// has to do so.
        /// </remarks>
        public DateTime WriteTime {
            get {
                return Utils.GetTimestepFileWriteTime(this);
            }
            set {
                realTimestepInfo.Value.WriteTime = value;
            }
        }

        /// <summary>
        /// See <see cref="TimestepInfo.Name"/>
        /// </summary>
        public string Name {
            get {
                return realTimestepInfo.Value.Name;
            }
            set {
                realTimestepInfo.Value.Name = value;
            }
        }

        /// <summary>
        /// The associated database
        /// </summary>
        public IDatabaseInfo Database {
            get {
                return Session.Database;
            }
        }

        /// <summary>
        /// See <see cref="TimestepInfo.CopyFor"/>
        /// </summary>
        /// <param name="targetDatabase"></param>
        /// <returns></returns>
        public ITimestepInfo CopyFor(IDatabaseInfo targetDatabase) {
            return realTimestepInfo.Value.CopyFor(targetDatabase);
        }

        #endregion
    }
}
