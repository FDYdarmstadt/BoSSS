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
using System.Runtime.Serialization;
using BoSSS.Foundation.Grid;
using ilPSP;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// interface definition of the database driver
    /// </summary>
    public interface IDatabaseDriver : IDisposable {

        /// <summary>
        /// the file-system driver
        /// </summary>
        IFileSystemDriver FsDriver {
            get;
        }

        //IFormatter Formatter {
        //    get;
        //}

        /// <summary>
        /// MPI rank of actual process within the MPI world communicator
        /// </summary>
        int MyRank {
            get;
        }

        /// <summary>
        /// Number of MPI processes within the MPI world communicator
        /// </summary>
        int Size {
            get;
        }

        /// <summary>
        /// Creates a new session;
        /// </summary>
        SessionInfo CreateNewSession(IDatabaseInfo database);

        /// <summary>
        /// tracing setup
        /// </summary>
        void InitTraceFile(SessionInfo si);

        ///// <summary>
        ///// tracing setup.
        ///// </summary>
        //void CloseTraceFile();

        /// <summary>
        /// Returns a write-stream for some new log file.
        /// </summary>
        Stream GetNewLogStream(SessionInfo si, string name);      

        /// <summary>
        /// saves a vector to the database, allocating a guid automatically
        /// </summary>
        /// <param name="vector"></param>
        /// <returns>
        /// the guid that was allocated to identify the vector within the storage system
        /// </returns>
        Guid SaveVector<T>(IList<T> vector);

        /// <summary>
        /// saves a vector to the database, under a specified Guid
        /// </summary>
        /// <param name="vector">
        /// the part of the vector which is stored on the local process.
        /// </param>
        /// <param name="id">
        /// the Guid under which the vector should be stored; must be the same
        /// on all MPI processes.
        /// </param>
        void SaveVector<T>(IList<T> vector, Guid id);

        /// <summary>
        /// Loads a vector from the database
        /// </summary>
        /// <param name="id"></param>
        /// <param name="part">
        /// Optional partition of the vector among MPI processors: if null, a partition is defined by the loader logic.
        /// </param>
        IList<T> LoadVector<T>(Guid id, ref Partitioning part);

        /// <summary>
        /// tests whether a grid with GUID <paramref name="g"/> exists in database, or not;
        /// </summary>
        bool GridExists(Guid g);

        /// <summary>
        /// Creates the grid info object for the grid with
        /// <paramref name="gridId"/> within the given
        /// <paramref name="database"/>.
        /// </summary>
        /// <param name="gridId"></param>
        /// <param name="database"></param>
        /// <returns></returns>
        IGridInfo LoadGridInfo(Guid gridId, IDatabaseInfo database);

        /// <summary>
        /// loads the grid identified by <paramref name="gridId"/> from the
        /// given database.
        /// </summary>
        /// <param name="gridId">The unique identifier of the grid.</param>
        /// <param name="database">
        /// The database that is associated with the grid.
        /// </param>
        /// <returns>
        /// The loaded grid
        /// </returns>
        IGrid LoadGrid(Guid gridId, IDatabaseInfo database);

        /// <summary>
        /// Loads the actual grid data for the given <paramref name="grid"/>.
        /// That is, loads the actual cell data.
        /// </summary>
        /// <param name="grid"></param>
        /// <returns></returns>
        IGrid LoadGridData(Grid.Classic.GridCommons grid);

        /// <summary>
        /// Loads the given <paramref name="sessionId"/> from the given
        /// <paramref name="database"/>.
        /// </summary>
        /// <param name="sessionId"></param>
        /// <param name="database"></param>
        /// <returns></returns>
        SessionInfo LoadSession(Guid sessionId, IDatabaseInfo database);

        /// <summary>
        /// Saves a session info object to a file on the disk.
        /// </summary>
        /// <param name="session">The session to be saved.</param>
        void SaveSessionInfo(ISessionInfo session);

        /// <summary>
        /// Searches for an equivalent grid in the database and, if none is found
        /// saves a grid object to the database.
        /// </summary>
        /// <param name="grd">
        /// On entry, the grid which should be saved to the database.
        /// On exit, either unchanged, or the equivalent grid.
        /// </param>
        /// <param name="EquivalentGridFound">
        /// Inidicates that an equivalent grid was found.
        /// </param>
        /// <param name="database"></param>
        Guid SaveGridIfUnique(ref IGrid grd, out bool EquivalentGridFound, IDatabaseInfo database);

        /// <summary>
        /// saves the grid object to the database;
        /// </summary>
        /// <returns>
        /// the Guid of the grid
        /// </returns>
        /// <param name="grd">
        /// the grid to save
        /// </param>
        /// <param name="database"></param>
        Guid SaveGrid(IGrid grd, IDatabaseInfo database);

        /// <summary>
        /// Saves a time-step to the database's persistent memory.
        /// </summary>
        /// <param name="physTime">Physical time of the time-step.</param>
        /// <param name="TimestepNo">Time-step number.</param>
        /// <param name="currentSession">The session associated with the time-step.</param>
        /// <param name="fields">The fields of the time-step.</param>
        /// <param name="g">grid data object (required if <paramref name="fields"/> is empty)</param>
        /// <returns>An object containing information about the time-step.</returns>
        TimestepInfo SaveTimestep(double physTime, TimestepNumber TimestepNo,
            SessionInfo currentSession, IGridData g, IEnumerable<DGField> fields);

        /// <summary>
        /// loads a single <see cref="TimestepInfo"/>-object from the database.
        /// </summary>
        TimestepInfo LoadTimestepInfo(Guid timestepGuid, ISessionInfo session, IDatabaseInfo database);

        /// <summary>
        /// Gathers all time-step IDs of a session.
        /// </summary>
        /// <param name="sessionGuid">ID of the session.</param>
        /// <returns>A collection of th session's time-step IDs.</returns>
        IEnumerable<Guid> GetTimestepGuids(Guid sessionGuid);

        /// <summary>
        /// Removes the given <paramref name="timestepGuid"/> from the
        /// time-step log for the given <paramref name="sessionGuid"/>
        /// </summary>
        /// <param name="sessionGuid"></param>
        /// <param name="timestepGuid"></param>
        void RemoveTimestepGuid(Guid sessionGuid, Guid timestepGuid);

        /// <summary>
        /// Loads a time-step from the database into previously allocated DG-fields (<paramref name="PreAllocatedFields"/>).
        /// </summary>
        void LoadFieldData(ITimestepInfo info, IGridData grdDat, IEnumerable<DGField> PreAllocatedFields);

        /// <summary>
        /// Loads a time-step from the database.
        /// </summary>
        /// <remarks>
        /// By using this method, it is ensured that the loaded/returned fields have the same DG polynomial degree as in the database.
        /// </remarks>
        IEnumerable<DGField> LoadFields(ITimestepInfo info, IGridData grdDat, IEnumerable<string> NameFilter = null);


    }
}
