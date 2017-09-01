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

using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace ilPSP.Utils {

    /// <summary>
    /// Wrapper for <see cref="IEnumerable{T}"/> that simplifies access to
    /// some handy information about the individual entities, such as the
    /// current index. This is mainly sensible in for-each loops where you are
    /// too lazy to add an additional counter or the like.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <remarks>
    /// Inspired by http://www.yoda.arachsys.com/csharp/miscutil/usage/smartenumerable.html
    /// </remarks>
    public class SmartEnumerable<T> : IEnumerable<SmartEnumerable<T>.Entry> {

        /// <summary>
        /// The wrapped enumerable
        /// </summary>
        private readonly IEnumerable<T> entities;

        /// <summary>
        /// Wraps the given enumerable
        /// </summary>
        /// <param name="entities">
        /// Entities to be enumerated smartly
        /// </param>
        public SmartEnumerable(IEnumerable<T> entities) {
            this.entities = entities;
        }

        #region IEnumerable<T> Members

        /// <summary>
        /// Wraps each entity into an <see cref="Entry"/> and returns the
        /// result as an enumerable. Note that the semantics of this method are
        /// slightly different from most of the standard methods related to
        /// <see cref="IEnumerable{T}"/> since we count the items eagerly
        /// (i.e., we enumerate over the full list before returning the first
        /// item)
        /// </summary>
        /// <returns></returns>
        public IEnumerator<Entry> GetEnumerator() {
            int lastIndex = entities.Count() - 1;
            return entities.Select((value, i) =>
                new Entry(value, i, i == 0, i == lastIndex)).GetEnumerator();
        }

        #endregion

        #region IEnumerable Members

        /// <summary>
        /// See <see cref="SmartEnumerable{T}.GetEnumerator"/>
        /// </summary>
        /// <returns></returns>
        IEnumerator IEnumerable.GetEnumerator() {
            return this.GetEnumerator();
        }

        #endregion

        /// <summary>
        /// An individual entry to be enumerated
        /// </summary>
        public struct Entry {

            /// <summary>
            /// Wraps the given data
            /// </summary>
            /// <param name="value">
            /// The actual entry
            /// </param>
            /// <param name="index">
            /// The index of the entry within the given enumerable
            /// </param>
            /// <param name="isFirst">
            /// True, if the given entry is the first entry in the given
            /// enumerable
            /// </param>
            /// <param name="isLast">
            /// True, if the given entry is the last entry in the given
            /// enumerable
            /// </param>
            public Entry(T value, int index, bool isFirst, bool isLast) {
                this.Value = value;
                this.Index = index;
                this.IsFirst = isFirst;
                this.IsLast = isLast;
            }

            /// <summary>
            /// The actual entry
            /// </summary>
            public readonly T Value;

            /// <summary>
            /// The index of the entry within the given enumerable
            /// </summary>
            public readonly int Index;

            /// <summary>
            /// True, if the given entry is the first entry in the given
            /// enumerable
            /// </summary>
            public readonly bool IsFirst;

            /// <summary>
            /// True, if the given entry is the last entry in the given
            /// enumerable
            /// </summary>
            public readonly bool IsLast;
        }
    }
}
