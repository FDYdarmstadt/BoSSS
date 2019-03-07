using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MPI.Wrappers;
using System.IO;
using ilPSP.Tracing;
using System.Diagnostics;

namespace BoSSS.Foundation.IO
{
    class SessionDatabaseDriver : MPIProcess, IDisposable
    {
        readonly VectorDataSerializer Driver;
        public SessionDatabaseDriver(VectorDataSerializer driver)
        {
            Driver = driver;
        }

        public void Dispose()
        {
            Driver.Dispose();
            if (m_stdout != null)
            {
                Debug.Assert(ilPSP.Environment.StdOut.WriterS.Contains(m_stdout));
                Debug.Assert(ilPSP.Environment.StdErr.WriterS.Contains(m_stderr));

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
        public SessionInfo CreateNewSession(IDatabaseInfo database)
        {

            Guid id = Guid.NewGuid();
            if (Driver.FsDriver == null)
                id = Guid.Empty;
            id = id.MPIBroadcast(0);

            // init driver
            // ===========
            if (Driver.FsDriver != null)
            {

                if (MyRank == 0)
                    Driver.FsDriver.CreateSessionDirectory(id);

                // ensure that the session directory is available, before ANY mpi process continues.
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                System.Threading.Thread.Sleep(1000);
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            }

            SessionInfo si = new SessionInfo(id, database);

            // create copy of stdout and stderr
            // ================================
            if (this.Driver.FsDriver != null)
            {
                m_stdout = this.Driver.FsDriver.GetNewLog("stdout." + MyRank, id);
                m_stderr = this.Driver.FsDriver.GetNewLog("stderr." + MyRank, id);

                ilPSP.Environment.StdOut.WriterS.Add(m_stdout);
                ilPSP.Environment.StdErr.WriterS.Add(m_stderr);
            }
            return si;
        }

        /// <summary>
        /// Saves a session info object to a file on the disk.
        /// </summary>
        /// <param name="session">The session to be saved.</param>
        public void SaveSessionInfo(ISessionInfo session)
        {
            using (Stream s = Driver.FsDriver.GetSessionInfoStream(true, session.ID))
            {
                Driver.Serialize(s, session);
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
        public SessionInfo LoadSession(Guid sessionId, IDatabaseInfo database)
        {
            using (var tr = new FuncTrace())
            {
                tr.Info("Loading session " + sessionId);

                using (Stream s = Driver.FsDriver.GetSessionInfoStream(false, sessionId))
                {
                    SessionInfo loadedSession = Driver.Deserialize<SessionInfo>(s);
                    loadedSession.Database = database;
                    loadedSession.WriteTime = Utils.GetSessionFileWriteTime(loadedSession);
                    s.Close();

                    return loadedSession;
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
        public static string GetSessionDirectory(ISessionInfo session)
        {
            string path = Path.Combine(
                session.Database.Path,
                StandardFsDriver.SessionsDir,
                session.ID.ToString());
            return path;
        }
    }
}
