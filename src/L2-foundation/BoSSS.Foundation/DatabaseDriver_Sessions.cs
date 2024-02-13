﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MPI.Wrappers;
using System.IO;
using ilPSP.Tracing;
using System.Diagnostics;

namespace BoSSS.Foundation.IO {
    class SessionDatabaseDriver : MPIProcess, IDisposable {
        readonly ISerializer Driver;
        IFileSystemDriver fsDriver;
        public SessionDatabaseDriver(ISerializer driver, IFileSystemDriver FsDriver) {
            fsDriver = FsDriver;
            Driver = driver;
        }

        public void Dispose() {
            if(m_stdout != null) {
                //Debug.Assert(ilPSP.Environment.StdOut.WriterS.Contains(m_stdout));
                //Debug.Assert(ilPSP.Environment.StdErr.WriterS.Contains(m_stderr));

                Console.Out.Flush();
                Console.Error.Flush();

                ilPSP.Environment.StdOut.WriterS.Remove(m_stdout);
                ilPSP.Environment.StdErr.WriterS.Remove(m_stderr);

                m_stderr.Close();
                m_stdout.Close();
                m_stderr.Dispose();
                m_stdout.Dispose();
            }
        }

        TextWriter m_stdout = null;
        TextWriter m_stderr = null;

        /// <summary>
        /// Creates a new session;
        /// </summary>
        public SessionInfo CreateNewSession(IDatabaseInfo database) {

            Guid id = Guid.NewGuid();
            if(fsDriver == null)
                id = Guid.Empty;
            id = id.MPIBroadcast(0);

            // init driver
            // ===========
            if(fsDriver != null) {

                if(MyRank == 0)
                    fsDriver.CreateSessionDirectory(id);

                // ensure that the session directory is available, before ANY mpi process continues.
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                System.Threading.Thread.Sleep(1000);
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            }

            SessionInfo si = new SessionInfo(id, database);

            // create copy of stdout and stderr
            // ================================
            if(this.fsDriver != null) {
                m_stdout = this.fsDriver.GetNewLog("stdout." + MyRank, id);
                m_stderr = this.fsDriver.GetNewLog("stderr." + MyRank, id);

                ilPSP.Environment.StdOut.WriterS.Add(m_stdout);
                ilPSP.Environment.StdErr.WriterS.Add(m_stderr);
            }
            return si;
        }

        /// <summary>
        /// Saves a session info object to a file on the disk.
        /// </summary>
        /// <param name="session">The session to be saved.</param>
        public void SaveSessionInfo(ISessionInfo session) {
            using(Stream s = fsDriver.GetSessionInfoStream(true, session.ID)) {
                Driver.Serialize(s, session, typeof(SessionInfo));
                s.Close();
            }
        }

        /// <summary>
        /// Loads the given <paramref name="sessionId"/> from the given
        /// <paramref name="database"/>.
        /// </summary>
        /// <param name="sessionId"></param>
        /// <param name="database"></param>
        /// <returns></returns>
        public SessionInfo LoadSession(Guid sessionId, IDatabaseInfo database) {
            using(var tr = new FuncTrace()) {
                tr.Info("Loading session " + sessionId);

                using(Stream s = fsDriver.GetSessionInfoStream(false, sessionId)) {
                    // fk, 19aug23: Create a MemoryStream and write he buffer into it
                    // for certain streams, the JSON deserialization ifrom a memory stream than much faster than directly deserializing from the original stream
                    using (MemoryStream memoryStream = new MemoryStream()) {
                        s.CopyTo(memoryStream);
                        memoryStream.Seek(0, SeekOrigin.Begin);
                        s.Close();

                        SessionInfo loadedSession = (SessionInfo)Driver.Deserialize(memoryStream, typeof(SessionInfo));
                        loadedSession.Database = database;
                        loadedSession.WriteTime = Utils.GetSessionFileWriteTime(loadedSession);

                        return loadedSession;
                    }
                }
            }
        }

        /// <summary>
        /// Retrieves the directory where the files for the selected
        /// <paramref name="session"/> are stored.
        /// </summary>
        /// <param name="session">
        /// The selected session.
        /// </param>
        /// <remarks>
        /// Should work on any System.
        /// </remarks>
        public static string GetSessionDirectory(ISessionInfo session) {
            string path = Path.Combine(
                session.Database.Path,
                StandardFsDriver.SessionsDir,
                session.ID.ToString());
            return path;
        }
    }
}
