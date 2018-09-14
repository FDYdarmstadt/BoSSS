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
    /// Manages operations on higher-level objects of a database,
    /// e.g. session and grid info objects, as opposed to the file-based
    /// operations in <see cref="BoSSS.Foundation.IO.IFileSystemDriver "/>
    /// </summary>
    public interface IDatabaseController {

        /// <summary>
        /// The driver's database
        /// </summary>
        IDatabaseInfo Database {
            get;
        }

        /// <summary>
        /// IO DatabaseDriver for lower-level operations
        /// </summary>
        IDatabaseDriver DBDriver {
            get;
        }

        /// <summary>
        /// The sessions associated with this driver's database
        /// </summary>
        IEnumerable<ISessionInfo> Sessions {
            get;
        }

        /// <summary>
        /// The grids associated with this driver's database
        /// </summary>
        IEnumerable<IGridInfo> Grids {
            get;
        }

        /// <summary>
        /// Obtains the session information associated with the given
        /// <paramref name="sessionID"/> through deserialization.
        /// </summary>
        /// <param name="sessionID">
        /// The id of the requested session
        /// </param>
        /// <returns></returns>
        ISessionInfo GetSessionInfo(Guid sessionID);

        /// <summary>
        /// Retrieves all sessions using a particular <paramref name="grid"/>.
        /// </summary>
        /// <param name="grid">The grid in question.</param>
        /// <returns>A collection of sessions using the grid.</returns>
        IEnumerable<ISessionInfo> GetSessionInfos(IGridInfo grid);

        /// <summary>
        /// Saves a session info object to the persistent memory.
        /// </summary>
        /// <param name="session">The session to be saved.</param>
        void SaveSessionInfo(ISessionInfo session);

        /// <summary>
        /// Deletes a session from a database.
        /// </summary>
        /// <param name="session">
        /// The session to be deleted.
        /// </param>
        void DeleteSession(ISessionInfo session);

        /// <summary>
        /// Copies a session to another database.
        /// </summary>
        /// <param name="session">The session to be copied.</param>
        /// <param name="dest">The destination database.</param>
        /// <returns>An info object for the session in the destination database.</returns>
        ISessionInfo CopySession(ISessionInfo session, IDatabaseInfo dest);

        /// <summary>
        /// Moves a session to another database.
        /// </summary>
        /// <param name="session">The session to be moved.</param>
        /// <param name="dest">The destination database.</param>
        /// <returns>An info object for the session in the destination database.</returns>
        ISessionInfo MoveSession(ISessionInfo session, IDatabaseInfo dest);
        
        /// <summary>
        /// Deletes a time-step from the database and its files from the file
        /// system.
        /// </summary>
        /// <param name="timestep">The time-step to be deleted.</param>
        void DeleteTimestep(ITimestepInfo timestep);

        /// <summary>
        /// Gets all time step info objects for a given session.
        /// </summary>
        /// <param name="session">
        /// The session in question.
        /// </param>
        /// <returns>
        /// A list of time-steps associated with <paramref name="session"/>.
        /// </returns>
        IList<ITimestepInfo> GetTimestepInfos(ISessionInfo session);

        /// <summary>
        /// Retrieves the grid with the given <paramref name="gridID"/>.
        /// </summary>
        /// <param name="gridID">The grid id.</param>
        /// <returns>
        /// The grid info object for the selected grid.
        /// </returns>
        IGridInfo GetGridInfo(Guid gridID);

        /// <summary>
        /// Adds an initialization context for a particular grid object, i.e.
        /// such the grid data object and grid-global objects (like a basis)
        /// don't have to be created over and over again.
        /// </summary>
        /// <param name="gridData"></param>
        void AddGridInitializationContext(IGridData gridData);

        /// <summary>
        /// Retrieves the initialization for a particular time-step.
        /// </summary>
        /// <param name="ts"></param>
        /// <returns></returns>
        IInitializationContext GetInitializationContext(ITimestepInfo ts);

        /// <summary>
        /// Saves a grid info object to the persistent memory.
        /// </summary>
        /// <param name="grid">The grid to be saved.</param>
        void SaveGridInfo(IGridInfo grid);

        /// <summary>
        /// Retrieves information about all grids used by a given session.
        /// </summary>
        /// <param name="session">
        /// The session in question.
        /// </param>
        /// <returns>
        /// The collection of grids associated with <paramref name="session"/>.
        /// </returns>
        IEnumerable<IGridInfo> GetGridInfos(ISessionInfo session);

        /// <summary>
        /// Retrieves all files associated with a grid.
        /// </summary>
        /// <param name="grid">The info object of the grid</param>
        /// <returns>Paths to all files associated with <paramref name="grid"/></returns>
        IEnumerable<string> GetGridFiles(IGridInfo grid);

        /// <summary>
        /// Deletes a grid from the database,
        /// </summary>
        /// <param name="grid">The grid to be deleted.</param>
        /// <param name="safelyDelete">
        /// true if should only be deleted if it is not in use by any session;
        /// false if it should be deleted regardless of its usage status.
        /// </param>
        /// <returns>true if the deletion was successful, false otherwise.</returns>
        bool DeleteGrid(IGridInfo grid, bool safelyDelete = true);

        /// <summary>
        /// Copies a grid to the given <paramref name="destination"/>
        /// </summary>
        /// <param name="grid">
        /// The grid to be copied
        /// </param>
        /// <param name="destination">
        /// The destination database
        /// </param>
        /// <returns>
        /// The grid info object for the copy of <paramref name="grid"/>.
        /// </returns>
        IGridInfo CopyGrid(IGridInfo grid, IDatabaseInfo destination);

        /// <summary>
        /// Clears the entire database, leaving only the basic folder
        /// structure behind
        /// </summary>
        void ClearDatabase();

        /// <summary>
        /// Disposes of all unused objects in the database
        /// </summary>
        void CleanDatabase();
    }
}