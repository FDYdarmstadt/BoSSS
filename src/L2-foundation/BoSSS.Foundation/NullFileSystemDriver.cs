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
    /// Null object for <see cref="IFileSystemDriver"/>
    /// </summary>
    public class NullFileSystemDriver : IFileSystemDriver {

        /// <summary>
        /// The one and only instance (reference equality)
        /// </summary>
        public static readonly NullFileSystemDriver Instance = new NullFileSystemDriver();

        /// <summary>
        /// Hides the public default constructor
        /// </summary>
        private NullFileSystemDriver() {
        }

        #region IFilesystemDriver Members

        /// <summary>
        /// Returns false.
        /// </summary>
        /// <param name="_SessionGuid"></param>
        /// <returns></returns>
        public bool IsSessionInUse(Guid _SessionGuid) {
            return false;
        }

        /// <summary>
        /// Returns the <see cref="StreamWriter.BaseStream"/> of
        /// <see cref="StreamWriter.Null"/>
        /// </summary>
        /// <param name="logName"></param>
        /// <param name="sessionGuid"></param>
        /// <returns></returns>
        public Stream GetNewLogStream(string logName, Guid sessionGuid) {
            return StreamWriter.Null.BaseStream;
        }

        /// <summary>
        /// Returns <see cref="TextWriter.Null"/>
        /// </summary>
        /// <param name="logName"></param>
        /// <param name="sessionGuid"></param>
        /// <returns></returns>
        public TextWriter GetNewLog(string logName, Guid sessionGuid) {
            return TextWriter.Null;
        }

        /// <summary>
        /// Returns an empty string
        /// </summary>
        /// <param name="sessionGuid"></param>
        /// <returns></returns>
        public string GetTimestepLogPath(Guid sessionGuid) {
            return string.Empty;
        }

        /// <summary>
        /// Not implemented.
        /// </summary>
        /// <param name="sessionGuid"></param>
        /// <returns></returns>
        public Stream GetTimestepLogStream(Guid sessionGuid) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Does nothing
        /// </summary>
        /// <param name="_SessionGuid"></param>
        public void CreateSessionDirectory(Guid _SessionGuid) {
        }

        /// <summary>
        /// Not implemented.
        /// </summary>
        /// <param name="create"></param>
        /// <param name="gridGuid"></param>
        /// <returns></returns>
        public Stream GetGridStream(bool create, Guid gridGuid) {
            throw new NotSupportedException(
                "Cannot load grid because no database was selected.");
        }

        /// <summary>
        /// Not implemented.
        /// </summary>
        /// <param name="create"></param>
        /// <param name="vecGuid"></param>
        /// <param name="part"></param>
        /// <returns></returns>
        public Stream GetDistVectorDataStream(bool create, Guid vecGuid, int part) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Not implemented.
        /// </summary>
        /// <param name="create"></param>
        /// <param name="id"></param>
        /// <returns></returns>
        public Stream GetTimestepStream(bool create, Guid id) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Not implemented.
        /// </summary>
        /// <param name="create"></param>
        /// <param name="sessionGuid"></param>
        /// <returns></returns>
        public Stream GetSessionInfoStream(bool create, Guid sessionGuid) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Returns an empty collection.
        /// </summary>
        /// <returns></returns>
        public IEnumerable<Guid> GetAllGridGUIDs() {
            return new List<Guid>();
        }

        /// <summary>
        /// Returns an empty list.
        /// </summary>
        /// <returns></returns>
        public IEnumerable<Guid> GetAllSessionGUIDs() {
            yield break;
        }

        /// <summary>
        /// Returns an empty list.
        /// </summary>
        /// <returns></returns>
        public IEnumerable<Guid> GetAllDataVectorGUIDs() {
            yield break;
        }

        /// <summary>
        /// Empty
        /// </summary>
        public string BasePath {
            get {
                return "";
            }
        }

        /// <summary>
        /// Always returns false
        /// </summary>
        public bool IsDisposed {
            get {
                return false;
            }
        }

        #endregion

        #region IDisposable Members

        /// <summary>
        /// Does nothing
        /// </summary>
        public void Dispose() {
        }

        #endregion
    }
}
