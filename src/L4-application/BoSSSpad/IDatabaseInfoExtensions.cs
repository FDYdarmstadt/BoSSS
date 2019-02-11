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

using BoSSS.Application.BoSSSpad;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.GridImport;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Extension methods for <see cref="IDatabaseInfo"/>
    /// </summary>
    public static class IDatabaseInfoExtensions {

        /// <summary>
        /// Deletes all files in this database.
        /// </summary>
        /// <param name="database">The database to be cleared.</param>
        public static void Clear(this IDatabaseInfo database) {
            Console.WriteLine("Database: " + database.ToString());
            Console.Write("This will delete _every_ file in the database. Do you wish to proceed? [y/n]: ");
            string line = Console.ReadLine();

            if (line.ToLower().Equals("y")) {
                database.Controller.ClearDatabase();
                Console.WriteLine("Database successfully cleared.");
            } else {
                Console.WriteLine("Database clear cancelled.");
            }
        }

        /// <summary>
        /// Deletes all unused files in this database.
        /// </summary>
        /// <param name="database">The database to be cleaned.</param>
        public static void Clean(this IDatabaseInfo database) {
            database.Controller.CleanDatabase();
        }

        /// <summary>
        /// Opens the database base directory.
        /// </summary>
        /// <param name="database">
        /// The selected database.
        /// </param>
        /// <remarks>
        /// Obviously, this only works in Windows environments.
        /// </remarks>
        public static void OpenBaseDirectory(this IDatabaseInfo database) {
            Process.Start(database.Path);
        }

        /// <summary>
        /// Imports the grid stored at the given location into
        /// <paramref name="database"/>.
        /// </summary>
        /// <param name="database">
        /// The database where the grid should be stored
        /// </param>
        /// <param name="fileName">
        /// The path to the grid file
        /// </param>
        /// <param name="test">
        /// If set to true, a <see cref="GridData"/> object is created in order to do more extensive testing on the imported grid.
        /// </param>
        /// <returns>
        /// The newly imported grid
        /// </returns>
        public static IGridInfo ImportGrid(this IDatabaseInfo database, string fileName, bool test = true) {
            GridCommons grid = GridImporter.Import(fileName);
            if (test)
                new GridData(grid);
            database.Controller.SaveGridInfo(grid);
            return grid;
        }

        /// <summary>
        /// Stores a grid in the database.
        /// </summary>
        /// <param name="database"></param>
        /// <param name="grd"></param>
        /// <returns></returns>
        public static Guid SaveGrid<TG>(this IDatabaseInfo database, ref TG grd) where TG : IGrid //
        {
            bool found;
            IGrid grid = grd;

            Guid GridGuid = database.Controller.DBDriver.SaveGridIfUnique(ref grid, out found, database);
            if (found) {
                Console.WriteLine("An equivalent grid is already present in the database -- the grid will not be saved.");
                grd = ((TG)grid);
            }
            return GridGuid;
        }


        /// <summary>
        /// Stores a grid in the database.
        /// </summary>
        /// <param name="database"></param>
        /// <param name="grd"></param>
        /// <returns></returns>
        public static IGrid SaveGrid(this IDatabaseInfo database, IGrid grd) {
            bool found;
            IGrid grid = grd;

            Guid GridGuid = database.Controller.DBDriver.SaveGridIfUnique(ref grid, out found, database);
            if (found) {
                Console.WriteLine("An equivalent grid is already present in the database -- the grid will not be saved.");
                grd = (GridCommons)grid;
            }
            return grid;
        }

        /// <summary>
        /// Save a list of DG fields (all defined on the same grid) into the database.
        /// </summary>
        /// <param name="database"></param>
        /// <param name="fields"></param>
        /// <returns></returns>
        public static ITimestepInfo SaveTimestep(this IDatabaseInfo database, params DGField[] fields) {
            return SaveTimestep<DGField>(database, fields);
        }

        /// <summary>
        /// Save a list of DG fields (all defined on the same grid) into the database.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="database"></param>
        /// <param name="fields"></param>
        /// <returns></returns>
        public static ITimestepInfo SaveTimestep<T>(this IDatabaseInfo database, IEnumerable<T> fields) where T : DGField {
            if (fields.Count() < 1)
                throw new ArgumentException("empty list");
            DGField[] _fields = fields.Select(f => (DGField)f).ToArray<DGField>();

            IGridData gDat = _fields[0].GridDat;
            for(int i = 0; i < _fields.Length; i++) {
                if (!object.ReferenceEquals(gDat, _fields[i].GridDat))
                    throw new ArgumentException("Mismatch of grids: DG field #" + i + " is defined on a different grid than the previous ones.");
            }

            var grid = gDat.Grid;
            if(grid.ID.Equals(Guid.Empty) || database.Grids.Where(g => g.ID.Equals(grid.ID)).Count() < 1) {
                // its necessary to save the grid to the database first

                /*
                var gridNew = SaveGrid(database, grid);
                if(!object.ReferenceEquals(gridNew, grid)) {
                    for (int i = 0; i < _fields.Length; i++) {
                        DGField fOld = _fields[i];
                        if(fOld.GetType() != typeof(SinglePhaseField)) {
                            throw
                        }
                    }
                }
                */

                // note: calling save grid would be easy,
                // but that would require re-creation of the DG fields on the new grid,
                // if 'SaveGrid' returns some *other* grid object that is already present in the database.

                throw new ArgumentException("Grid is not stored in database; Save the grid to database in order to use this method.");
            }

            SessionInfo dummySession = new SessionInfo(Guid.NewGuid(), database);
            dummySession.Name = "InitialValueSession";
            dummySession.ProjectName = InteractiveShell.WorkflowMgm.CurrentProject;
            dummySession.Save();
            return database.Controller.DBDriver.SaveTimestep(0.0, 0, dummySession, gDat , _fields);
        }
    }
}
