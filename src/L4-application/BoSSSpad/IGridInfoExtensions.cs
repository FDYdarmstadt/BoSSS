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

using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using BoSSS.Application.BoSSSpad;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;


namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Extensions methods for <see cref="IGridInfo"/>.
    /// </summary>
    public static class IGridInfoExtensions {

        /// <summary>
        /// Copies a session to another database.
        /// </summary>
        /// <param name="grid">The grid to be copied.</param>
        /// <param name="targetDB">The target database.</param>
        public static void Copy(this IGridInfo grid, IDatabaseInfo targetDB) {
            grid.Database.Controller.CopyGrid(grid, targetDB);
        }

        /// <summary>
        /// Deletes a grid from its database.
        /// </summary>
        /// <param name="grid">The grid to be copied.</param>
        /// <param name="force">
        /// If false, the grid is only deleted if it is not in use.
        /// If true, the deletion will be forced.
        /// </param>
        public static void Delete(this IGridInfo grid, bool force = false) {
            // Only ask for confirmation of unsafe deletion.
            if (force == false) {
                Console.WriteLine("Grid: " + grid.ToString());
                Console.Write("The grid will be removed, even if it is in use by a session. Do you wish to proceed? [y/n]: ");
                string line = Console.ReadLine();

                if (line.ToLower().Equals("y")) {
                    grid.Database.Controller.DeleteGrid(grid, false);
                    Console.WriteLine("Grid " + grid.ID.ToString() + " deleted.");
                } else {
                    Console.WriteLine("Grid delete canceled.");
                }
            } else {
                // Don't need to ask for confirmation of safe delete.
                bool outcome = grid.Database.Controller.DeleteGrid(grid, true);
                if (outcome == true) {
                    Console.WriteLine("Grid " + grid.ID.ToString() + " deleted.");
                } else {
                    Console.WriteLine("Grid was not deleted. Probably still in use by a session.");
                }
            }
        }

        /// <summary>
        /// Deletes all grids inside a collection, provided that they can be
        /// deleted. The user is not asked to confirm this action once, and
        /// not for every individual element in the collection.
        /// </summary>
        /// <param name="grids">The entities to be deleted.</param>
        public static void DeleteAll(this IEnumerable<IGridInfo> grids) {
            if (grids.IsNullOrEmpty()) {
                Console.WriteLine("Given collection is empty; nothing to delete");
                return;
            }

            Console.WriteLine("Grids to delete:");
            foreach (IGridInfo grid in grids) {
                Console.WriteLine(grid.ToString());
            }
            Console.Write("Do you really want to delete these grids? [y/n]: ");
            string line = Console.ReadLine();

            if (line.ToLower().Equals("y")) {
                foreach (IGridInfo grid in grids) {
                    grid.Database.Controller.DeleteGrid(grid);
                }
                Console.WriteLine("Grids successfully deleted.");
            } else {
                Console.WriteLine("Grid delete canceled.");
            }
        }

        /// <summary>
        /// Experimental
        /// </summary>
        /// <param name="grid"></param>
        public static void Save(this IGridInfo grid) {
            if (grid is GridCommons) {
                GridCommons realGrid = (GridCommons)grid;
                grid.Database.Controller.DBDriver.SaveGrid(realGrid);
            } else {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Returns the single grid that matches the given
        /// <paramref name="guid"/>.
        /// </summary>
        /// <param name="grids">
        /// The list of grids to be queried.
        /// </param>
        /// <param name="guid">
        /// The guid of the grid in question.
        /// </param>
        /// <returns>
        /// A single grid that matches <paramref name="guid"/>, if it
        /// exists. Otherwise, an error will be thrown.
        /// </returns>
        public static IGridInfo Find(this IEnumerable<IGridInfo> grids, PartialGuid guid) {
            return grids.Where(g => g.ID == guid).Single();
        }

        /// <summary>
        /// Opens the directory where the export for the selected
        /// <paramref name="grid"/> are stored in the explorer.
        /// </summary>
        /// <param name="grid">
        /// The selected grid.
        /// </param>
        /// <remarks>
        /// Obviously, this only works in Windows environments.
        /// </remarks>
        public static void OpenExportDirectory(this IGridInfo grid) {
            Process.Start(Utils.GetExportDirectory(grid));
        }

        /// <summary>
        /// Prints the directory where the exports for the selected
        /// <paramref name="grid"/> are stored to the console.
        /// </summary>
        /// <param name="grid">
        /// The selected grid.
        /// </param>
        /// <remarks>
        /// Should work on any System.
        /// </remarks>
        public static void PrintExportDirectory(this IGridInfo grid) {
            Console.WriteLine(Utils.GetExportDirectory(grid));
        }

        /// <summary>
        /// Convenience interface to create a
        /// <see cref="GridExportInstruction"/>
        /// </summary>
        /// <param name="grid"></param>
        /// <returns></returns>
        public static GridExportInstruction Export(this IGridInfo grid) {
            return new GridExportInstruction(grid);
        }

        /// <summary>
        /// Retrieves all sessions using this grid.
        /// </summary>
        /// <param name="grid">The grid in question.</param>
        /// <returns>A collection of sessions associated with the grid.</returns>
        public static IEnumerable<ISessionInfo> GetSessions(this IGridInfo grid) {
            return grid.Database.Controller.GetSessionInfos(grid);
        }

        /// <summary>
        /// Determines the mesh size of a given <paramref name="gridInfo"/>
        /// based on <see cref="GridData.CellData.h_maxGlobal"/>.
        /// </summary>
        /// <param name="gridInfo">A grid</param>
        /// <returns>
        /// A (heuristic) measure for the mesh size of
        /// <paramref name="gridInfo"/>.
        /// </returns>
        public static double GetMeshSize(this IGridInfo gridInfo) {
            if (gridInfo is GridProxy) {
                gridInfo = gridInfo.Cast<GridProxy>().RealGrid;
            }

            GridCommons grid = gridInfo as GridCommons;
            if (!(grid is GridCommons)) {
                throw new ArgumentException("Only works for grids of type 'GridCommons'");
            }

            return new GridData(grid).Cells.h_maxGlobal;
        }
    }
}
