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
using System.Linq;
using Renci.SshNet;

namespace BoSSS.Foundation.IO {

    public class SFTPFileSystemDriver : IFileSystemDriver, IDisposable {

        private SftpClient sftp;

        private string path;

        private bool isInitialized = false;

        public SFTPFileSystemDriver(ConnectionInfo connectionInfo, string path) {
            sftp = new SftpClient(connectionInfo);

            this.ConnectionInfo = connectionInfo;
            this.path = path;
        }

        public ConnectionInfo ConnectionInfo {
            get;
            private set;
        }

        private void EnsureIsInitialized() {
            if (IsDisposed) {
                throw new InvalidOperationException();
            }

            if (!isInitialized) {
                sftp.Connect();

                if (!sftp.Exists(path)) {
                    throw new ArgumentException(
                        "Database Error: base directory (" + path + ") does not exist.");
                }

                VerifyDirectoryStructure();
                isInitialized = true;
            }
        }

        #region IFileSystemDriver Members

        public Stream GetNewLogStream(string logName, Guid sessionGuid) {
            EnsureIsInitialized();
            string relPath = CombinePaths(
                StandardFsDriver.SessionsDir, sessionGuid.ToString(), logName + ".txt");
            return OpenFile(relPath, FileAccess.ReadWrite, true);
        }

        public TextWriter GetNewLog(string logName, Guid sessionGuid) {
            EnsureIsInitialized();
            return new StreamWriter(GetNewLogStream(logName, sessionGuid));
        }

        public string GetTimestepLogPath(Guid sessionGuid) {
            EnsureIsInitialized();
            return CombinePaths(
                path,
                StandardFsDriver.SessionsDir,
                sessionGuid.ToString(),
                "TimestepLog.txt");
        }

        public Stream GetTimestepLogStream(Guid sessionGuid) {
            EnsureIsInitialized();
            string filename = GetTimestepLogPath(sessionGuid);
            return sftp.Open(filename, FileMode.Open, FileAccess.Read);
        }

        public void CreateSessionDirectory(Guid sessionGuid) {
            if (sessionGuid.Equals(Guid.Empty))
                throw new ArgumentException();

            EnsureIsInitialized();

            // create session directory
            string SessionDir = CombinePaths(
                path, StandardFsDriver.SessionsDir, sessionGuid.ToString());
            sftp.CreateDirectory(SessionDir);
        }

        public Stream GetGridStream(bool createIfNotExists, Guid gridGuid) {
            EnsureIsInitialized();
            string relPath = CombinePaths(StandardFsDriver.GridsDir, gridGuid + ".grid");
            return OpenFile(relPath, FileAccess.Read, false);
        }

        public Stream GetDistVectorDataStream(bool create, Guid vecGuid, int part) {
            EnsureIsInitialized();
            string relPath = CombinePaths(
                StandardFsDriver.DistVectorDataDir,
                vecGuid.ToString() + "." + (part + 1) + ".data");
            return OpenFile(relPath, FileAccess.Read, true);
        }

        public Stream GetTimestepStream(bool create, Guid id) {
            EnsureIsInitialized();
            string relPath = CombinePaths(
                StandardFsDriver.TimestepDir, id.ToString() + ".ts");
            return OpenFile(relPath, FileAccess.Read, create);
        }

        public Stream GetSessionInfoStream(bool create, Guid sessionGuid) {
            EnsureIsInitialized();
            string relPath = CombinePaths(
                StandardFsDriver.SessionsDir,
                sessionGuid.ToString(),
                "Session.info");

            return OpenFile(relPath, FileAccess.Read, create);
        }

        public IEnumerable<Guid> GetAllGridGUIDs() {
            EnsureIsInitialized();
            return sftp.ListDirectory(CombinePaths(path, StandardFsDriver.GridsDir)).
                Where(f => f.IsRegularFile).
                Where(f => f.Name.EndsWith(".grid")).
                Select(f => new Guid(f.Name.Substring(0, 36))).ToArray();
        }

        public IEnumerable<Guid> GetAllSessionGUIDs() {
            EnsureIsInitialized();
            return sftp.ListDirectory(CombinePaths(path, StandardFsDriver.SessionsDir)).
                Where(f => f.IsDirectory).
                Where(f => f.Length == 36).
                Select(f => new Guid(f.Name)).ToArray();
        }

        public IEnumerable<Guid> GetAllDataVectorGUIDs() {
            EnsureIsInitialized();
            return sftp.ListDirectory(CombinePaths(path, StandardFsDriver.DistVectorDataDir)).
                Where(f => f.IsRegularFile).
                Where(f => f.Name.EndsWith(".1.data")).
                Select(f => new Guid(f.Name.Substring(0, 36))).ToArray();
        }
        
        public string BasePath {
            get;
            private set;
        }

        public bool IsDisposed {
            get;
            private set;
        }

        #endregion

        #region IDisposable Members

        public void Dispose() {
            if (sftp != null) {
                sftp.Dispose();
                sftp = null;
                isInitialized = false;
            }
            IsDisposed = true;
        }

        #endregion

        private string CombinePaths(params string[] paths) {
            return paths.Aggregate((s, t) => s + "/" + t);
        }

        private Stream CreateFile(string fullPath, bool force = false) {
            EnsureIsInitialized();
            if (force) {
                return sftp.Open(fullPath, FileMode.Create);
            } else {
                return sftp.Open(fullPath, FileMode.CreateNew);
            }
        }

        private Stream OpenFile(string relPath, FileAccess fileAccess, bool createIfNotExists = false) {
            EnsureIsInitialized();

            string fullPath = CombinePaths(path, relPath);
            if (sftp.Exists(fullPath)) {
                return sftp.Open(fullPath, FileMode.Open, fileAccess);
            } else {
                if (createIfNotExists) {
                    return CreateFile(fullPath, false);
                } else {
                    throw new FileNotFoundException(String.Format(
                        "File '{0}' does not exist", relPath));
                }
            }
        }

        private void VerifyDirectoryStructure() {
            List<string> directories = sftp.ListDirectory(path).
                Where(d => d.IsDirectory).
                Select(d => d.Name).
                ToList();

            if (!directories.Contains(StandardFsDriver.TimestepDir)) {
                throw new ArgumentException(
                    "Database Error: field state directory (" + StandardFsDriver.TimestepDir + ") does not exist.");
            }

            if (!directories.Contains(StandardFsDriver.DistVectorDataDir)) {
                throw new ArgumentException(
                    "Database Error: distributed vector data directory (" + StandardFsDriver.DistVectorDataDir + ") does not exist.");
            }

            if (!directories.Contains(StandardFsDriver.SessionsDir)) {
                throw new ArgumentException(
                    "Database Error: session directory (" + StandardFsDriver.SessionsDir + ") does not exist.");
            }

            if (!directories.Contains(StandardFsDriver.GridsDir)) {
                throw new ArgumentException(
                    "Database Error: grid directory (" + StandardFsDriver.GridsDir + ") does not exist.");
            }
        }
    }
}
