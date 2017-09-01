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
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// A proxy for <see cref="GridCommons"/> objects that allows for lazy
    /// loading
    /// </summary>
    public class GridProxy : IGridInfo {

        /// <summary>
        /// The real grid reflected by this object.
        /// </summary>
        private ExpirableLazy<Grid.Classic.GridCommons> realGrid;

        /// <summary>
        /// Indicates whether the grid data associated with
        /// <see cref="realGrid"/> has already been loaded. Note that this can
        /// be different from
        /// <see cref="ExpirableLazy{GridCommons}.IsValueCreated"/>
        /// </summary>
        private bool gridDataLoaded = false;

        /// <summary>
        /// Creates a proxy for the grid with id <paramref name="gridGuid"/>
        /// within the given <paramref name="database"/>
        /// </summary>
        /// <param name="gridGuid"></param>
        /// <param name="database"></param>
        public GridProxy(Guid gridGuid, IDatabaseInfo database) {
            this.ID = gridGuid;
            this.Database = database;
            realGrid = new ExpirableLazy<GridCommons>(
                () => database.Controller.DBDriver.LoadGridInfo(gridGuid, database).Cast<GridCommons>(),
                g => Utils.GetGridFileWriteTime(g) == g.WriteTime);
        }

        /// <summary>
        /// The real grid reflected by this object.
        /// </summary>
        /// <remarks>
        /// Induces a call to <see cref="EnsureGridDataIsLoaded"/>
        /// </remarks>
        public Grid.Classic.GridCommons RealGrid {
            get {
                EnsureGridDataIsLoaded();
                return realGrid.Value;
            }
        }

        /// <summary>
        /// See <see cref="GridCommons.ToString"/>
        /// </summary>
        /// <returns></returns>
        public override string ToString() {
            return realGrid.Value.ToString();
        }

        /// <summary>
        /// If <see cref="gridDataLoaded"/> is false, uses
        /// <see cref="IDatabaseDriver.LoadGridData"/> to load the actual cell
        /// data. Note that this is an expensive operation.
        /// </summary>
        private void EnsureGridDataIsLoaded() {
            if (!gridDataLoaded) {
                Database.Controller.DBDriver.LoadGridData(realGrid.Value);
            }
            gridDataLoaded = true;
        }

        #region IGridInfo Members

        /// <summary>
        /// See <see cref="GridCommons.NumberOfCells"/>
        /// </summary>
        public int NumberOfCells {
            get {
                return realGrid.Value.NumberOfCells;
            }
        }

        /// <summary>
        /// See <see cref="GridCommons.SpatialDimension"/>
        /// </summary>
        public int SpatialDimension {
            get {
                return realGrid.Value.SpatialDimension;
            }
        }

        #endregion

        #region IDatabaseEntityInfo<IGridInfo> Members

        /// <summary>
        /// The id of the grid.
        /// </summary>
        public Guid ID {
            get;
            private set;
        }

        /// <summary>
        /// See <see cref="GridCommons.CreationTime"/>
        /// </summary>
        public DateTime CreationTime {
            get {
                return realGrid.Value.CreationTime;
            }
        }

        /// <summary>
        /// See <see cref="GridCommons.WriteTime"/>.
        /// </summary>
        /// <remarks>
        /// Note that writing to this property induces a call to
        /// <see cref="EnsureGridDataIsLoaded"/>, while reading does not.
        /// </remarks>
        public DateTime WriteTime {
            get {
                return Utils.GetGridFileWriteTime(this);
            }
            set {
                EnsureGridDataIsLoaded();
                realGrid.Value.WriteTime = value;
            }
        }

        /// <summary>
        /// See <see cref="GridCommons.Name"/>.
        /// </summary>
        /// <remarks>
        /// Note that writing to this property induces a call to
        /// <see cref="EnsureGridDataIsLoaded"/>, while reading does not.
        /// </remarks>
        public string Name {
            get {
                return realGrid.Value.Name;
            }
            set {
                EnsureGridDataIsLoaded();
                realGrid.Value.Name = value;
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
        /// Creates a new grid proxy associated to database
        /// <paramref name="targetDatabase"/>.
        /// </summary>
        /// <param name="targetDatabase"></param>
        /// <returns></returns>
        public IGridInfo CopyFor(IDatabaseInfo targetDatabase) {
            return new GridProxy(ID, targetDatabase);
        }

        #endregion

        #region IEquatable<IGridInfo> Members

        /// <summary>
        /// See <see cref="GridCommons.Equals"/>
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals(IGridInfo other) {
            return realGrid.Value.Equals(other);
        }

        #endregion
    }
}
