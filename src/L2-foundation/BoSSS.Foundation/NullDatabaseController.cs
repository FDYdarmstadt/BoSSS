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

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// A null object for <see cref="IDatabaseController"/>
    /// </summary>
    public class NullDatabaseController : IDatabaseController {

        /// <summary>
        /// The one and only instance (reference equality)
        /// </summary>
        public static readonly NullDatabaseController Instance =
            new NullDatabaseController();

        /// <summary>
        /// Hides the public default constructor
        /// </summary>
        private NullDatabaseController() {
        }

        #region IDatabaseController Members

        /// <summary>
        /// Not implemented.
        /// </summary>
        public IDatabaseInfo Database {
            get {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Returns a new driver using a <see cref="NullFileSystemDriver"/>
        /// </summary>
        public IDatabaseDriver DBDriver {
            get {
                return new DatabaseDriver(NullFileSystemDriver.Instance);
            }
        }

        /// <summary>
        /// An empty list.
        /// </summary>
        public IEnumerable<ISessionInfo> Sessions {
            get {
                yield break;
            }
        }

        /// <summary>
        /// An empty list.
        /// </summary>
        public IEnumerable<IGridInfo> Grids {
            get {
                yield break;
            }
        }

        /// <summary>
        /// Returns null
        /// </summary>
        /// <param name="sessionID"></param>
        /// <returns></returns>
        public ISessionInfo GetSessionInfo(Guid sessionID) {
            return null;
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="gridData"></param>
        public void AddGridInitializationContext(IGridData gridData) {
        }

        /// <summary>
        /// Returns null.
        /// </summary>
        /// <param name="timestepGuid"></param>
        /// <returns></returns>
        public IInitializationContext GetInitializationContext(ITimestepInfo timestepGuid) {
            return null;
        }

        /// <summary>
        /// Returns an empty list
        /// </summary>
        /// <param name="grid"></param>
        /// <returns></returns>
        public IEnumerable<ISessionInfo> GetSessionInfos(IGridInfo grid) {
            yield break;
        }

        /// <summary>
        /// Does nothing
        /// </summary>
        /// <param name="session"></param>
        public void SaveSessionInfo(ISessionInfo session) {
        }

        /// <summary>
        /// Does nothing
        /// </summary>
        /// <param name="session"></param>
        public void DeleteSession(ISessionInfo session) {
        }

        /// <summary>
        /// Does nothing
        /// </summary>
        /// <param name="session"></param>
        /// <param name="dest"></param>
        /// <returns></returns>
        public ISessionInfo CopySession(ISessionInfo session, IDatabaseInfo dest) {
            return null;
        }

        /// <summary>
        /// Does nothing
        /// </summary>
        /// <param name="session"></param>
        /// <param name="dest"></param>
        /// <returns></returns>
        public ISessionInfo MoveSession(ISessionInfo session, IDatabaseInfo dest) {
            return null;
        }

        /// <summary>
        /// Does nothing
        /// </summary>
        /// <param name="timestep"></param>
        public void DeleteTimestep(ITimestepInfo timestep) {
        }

        /// <summary>
        /// Returns an empty list
        /// </summary>
        /// <param name="session"></param>
        /// <returns></returns>
        public IList<ITimestepInfo> GetTimestepInfos(ISessionInfo session) {
            return new List<ITimestepInfo>();
        }

        /// <summary>
        /// Returns null
        /// </summary>
        /// <param name="gridID"></param>
        /// <returns></returns>
        public IGridInfo GetGridInfo(Guid gridID) {
            return null;
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="grid"></param>
        public void SaveGridInfo(IGridInfo grid) {
        }

        /// <summary>
        /// Returns an empty list
        /// </summary>
        /// <param name="session"></param>
        /// <returns></returns>
        public IEnumerable<IGridInfo> GetGridInfos(ISessionInfo session) {
            yield break;
        }

        /// <summary>
        /// Returns an empty list
        /// </summary>
        /// <param name="grid"></param>
        /// <returns></returns>
        public IEnumerable<string> GetGridFiles(IGridInfo grid) {
            yield break;
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="grid"></param>
        /// <returns></returns>
        public bool DeleteGrid(IGridInfo grid) {
            return false;
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="grid"></param>
        /// <param name="safelyDelete"></param>
        /// <returns></returns>
        public bool DeleteGrid(IGridInfo grid, bool safelyDelete) {
            return false;
        }

        /// <summary>
        /// Does nothing and returns null.
        /// </summary>
        /// <param name="grid"></param>
        /// <param name="destination"></param>
        /// <returns></returns>
        public IGridInfo CopyGrid(IGridInfo grid, IDatabaseInfo destination) {
            return null;
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        public void ClearDatabase() {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        public void CleanDatabase() {
        }

        #endregion
    }
}
