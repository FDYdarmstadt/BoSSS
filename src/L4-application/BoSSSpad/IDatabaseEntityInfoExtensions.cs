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

using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Extension methods for <see cref="IDatabaseEntityInfo{T}"/>
    /// </summary>
    public static class IDatabaseEntityInfoExtensions {

        /// <summary>
        /// Retrieves the newest object from a list of
        /// <paramref name="entities"/>.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the individual entities
        /// </typeparam>
        /// <param name="entities">
        /// The entities in question
        /// </param>
        /// <returns>
        /// The entity that has been created last
        /// </returns>
        public static T Newest<T>(this IEnumerable<T> entities)
            where T : IDatabaseEntityInfo<T> {
            return entities.OrderByDescending(e => e.CreationTime).First();
        }

        /// <summary>
        /// Retrieves the <paramref name="count"/> newest objects from a list
        /// of <paramref name="entities"/>.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the individual entities
        /// </typeparam>
        /// <param name="entities">
        /// The entities in question
        /// </param>
        /// <param name="count">
        /// The number of entities to be retrieved
        /// </param>
        /// <returns>
        /// The <paramref name="count"/> entities within
        /// <paramref name="entities"/> that have been created last
        /// </returns>
        public static IEnumerable<T> Newest<T>(this IEnumerable<T> entities, int count)
            where T : IDatabaseEntityInfo<T> {
            return entities.OrderByDescending(e => e.CreationTime).Take(count);
        }

        /// <summary>
        /// Retrieves the oldest object from a list of
        /// <paramref name="entities"/>.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the individual entities
        /// </typeparam>
        /// <param name="entities">
        /// The entities in question
        /// </param>
        /// <returns>
        /// The entity that has been created first
        /// </returns>
        public static T Oldest<T>(this IEnumerable<T> entities)
            where T : IDatabaseEntityInfo<T> {
            return entities.OrderBy(e => e.CreationTime).First();
        }

        /// <summary>
        /// Retrieves the <paramref name="count"/> oldest object from a list of
        /// <paramref name="entities"/>.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the individual entities
        /// </typeparam>
        /// <param name="entities">
        /// The entities in question
        /// </param>
        /// <param name="count">
        /// The number of entities to be retrieved
        /// </param>
        /// <returns>
        /// The <paramref name="count"/> entities that have been created first
        /// </returns>
        public static IEnumerable<T> Oldest<T>(this IEnumerable<T> entities, int count)
            where T : IDatabaseEntityInfo<T> {
            return entities.OrderBy(e => e.CreationTime).Take(count);
        }

        /// <summary>
        /// Returns every <paramref name="nth"/> item in
        /// <paramref name="entities"/>.
        /// </summary>
        /// <param name="entities">The entities to be filtered</param>
        /// <param name="nth">The stride</param>
        /// <returns></returns>
        public static IEnumerable<T> Every<T>(this IEnumerable<T> entities, int nth)
            where T : IDatabaseEntityInfo<T> {
            return entities.Where((t, i) => i % nth == 0);
        }
    }
}
