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

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// This interface is used by the IO master (<see cref="DatabaseDriver"/>) to access some 
    /// storage device which supports streaming.
    /// </summary>
    public interface IFileSystemDriver : IDisposable {

        /// <summary>
        /// Returns a stream for logging information within this session;
        /// </summary>
        /// <param name="logName">
        /// name for the log file;
        /// this name identifies the log within the session;
        /// </param>
        /// <param name="sessionGuid">
        /// the session in which the log should be created.
        /// </param>        
        /// <returns></returns>
        Stream GetNewLogStream(string logName, Guid sessionGuid);

        /// <summary>
        /// returns a text writer for logging information within this session;
        /// </summary>
        /// <param name="logName">
        /// name for the log file;
        /// this name identifies the log within the session;
        /// </param>
        /// <param name="sessionGuid">
        /// the session in which the log should be created.
        /// </param>
        /// <returns></returns>
        TextWriter GetNewLog(string logName, Guid sessionGuid);

        /// <summary>
        /// Returns the path to the time-step log.
        /// </summary>
        /// <param name="sessionGuid">The ID of the session.</param>
        /// <returns>The path to the session's time-step log.</returns>
        string GetTimestepLogPath(Guid sessionGuid);

        /// <summary>
        /// Returns a text reader for reading the time-step log.
        /// </summary>
        /// <param name="sessionGuid">The ID of the session.</param>
        /// <returns>A text reader for the session's time-step log.</returns>
        Stream GetTimestepLogStream(Guid sessionGuid);

        /// <summary>
        /// creates a new session directory.
        /// </summary>
        void CreateSessionDirectory(Guid _SessionGuid);

        /// <summary>
        /// a stream to serialize grid objects (objects/classes derived from 
        /// <see cref="BoSSS.Foundation.Grid.GridCommons"/>);
        /// </summary>
        /// <param name="create">
        /// true for a write stream, false for a read stream
        /// </param>
        /// <param name="gridGuid">
        /// The id of the grid (see
        /// <see cref="BoSSS.Foundation.Grid.GridCommons.GridGuid"/>); If
        /// <paramref name="create"/> is true, the caller ensures that it is a
        /// newly created Guid. Otherwise, if <paramref name="create"/> is
        /// false and no field-state with specified Guid exists in the
        /// database, some exception should be thrown. The caller (some method
        /// of the IO-DatabaseDriver class <see cref="DatabaseDriver"/>)
        /// ensures that this parameter is equal on all MPI processes;
        /// </param>
        /// <returns></returns>
        Stream GetGridStream(bool create, Guid gridGuid);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="create">true for a write stream, false for a read stream</param>
        /// <param name="vecGuid">
        /// The caller (some method of the IO-DatabaseDriver class <see cref="DatabaseDriver"/>)
        /// ensures that this parameter is equal on all
        /// MPI processes;
        /// </param>
        /// <param name="part"></param>
        /// <returns></returns>
        Stream GetDistVectorDataStream(bool create, Guid vecGuid, int part);

        /// <summary>
        /// if no session guid was set (cf. <see cref="CreateSessionDirectory"/>), and
        /// <paramref name="create"/> is true, a call to this method should
        /// either return null or throw an exception.
        /// </summary>
        /// <param name="create">
        /// true for a write stream, false for a read stream
        /// </param>
        /// <param name="id">
        /// The id of the time step;
        /// If <paramref name="create"/> is true, the caller
        /// (<see cref="BoSSS.Foundation.IO.DatabaseDriver"/>) ensures that it
        /// is a newly created Guid. Otherwise, if <paramref name="create"/> is
        /// false and no time-step with specified Guid exists in the database,
        /// some exception should be thrown. The caller (some method of the
        /// IO-DatabaseDriver class <see cref="DatabaseDriver"/>) ensures that
        /// this parameter is equal on all MPI processes;
        /// </param>
        /// <returns></returns>
        Stream GetTimestepStream(bool create, Guid id);

        /// <summary>
        /// returns a stream to save/load "SessionInfo"-objects;
        /// </summary>
        /// <param name="create">
        /// true for a write stream, false for a read stream;
        /// </param>
        /// <param name="sessionGuid">
        /// the guid of the session for which the stream should be acquired;
        /// </param>
        /// <returns></returns>
        Stream GetSessionInfoStream(bool create, Guid sessionGuid);

        /// <summary>
        /// a list of all gird guids
        /// </summary>
        /// <returns></returns>
        IEnumerable<Guid> GetAllGridGUIDs();

        /// <summary>
        /// A list of the GUIDS of all sessions
        /// </summary>
        /// <returns></returns>
        IEnumerable<Guid> GetAllSessionGUIDs();

        /// <summary>
        /// a list of all Guid's of all data vectors;
        /// </summary>
        /// <returns></returns>
        IEnumerable<Guid> GetAllDataVectorGUIDs();

        /// <summary>
        /// Root path of this driver
        /// </summary>
        string BasePath {
            get;
        }

        /// <summary>
        /// Indicates whether this object has already been disposed
        /// </summary>
        bool IsDisposed {
            get;
        }
    }
}

