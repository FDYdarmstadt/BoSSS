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
using System.Linq;
using System.Text;
using ilPSP;
using BoSSS.Application.BoSSSpad;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Extension methods for <see cref="IEnumerable{T}"/>
    /// </summary>
    public static class IEnumerableExtensions {

        /// <summary>
        /// Format the given list of <paramref name="entities"/> for the
        /// console.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the entities.
        /// </typeparam>
        /// <param name="entities">
        /// The sequence of entities to be formatted.
        /// </param>
        /// <returns>
        /// A string representation of the list that is useful for console
        /// output.
        /// </returns>
        public static string Summary<T>(this IEnumerable<T> entities) {
            StringBuilder builder = new StringBuilder();

            builder.AppendLine(String.Format(
                "List of {0} instances of {1} {{",
                entities.Count(),
                typeof(T).GetGenericTypeName()));
            foreach (var entry in entities.AsSmartEnumerable()) {
                builder.AppendLine(String.Format(
                    "  {0}: {1}", entry.Index, entry.Value));
            }
            builder.Append("}");

            return builder.ToString();
        }

        /// <summary>
        /// Copies all entities inside a collection to another database,
        /// provided that they are database entities that can be copied,
        /// i.e. ISessionInfo and IGridInfo as of now.
        /// </summary>
        /// <typeparam name="T">The type of the entities.</typeparam>
        /// <param name="entities">The entities to be copied.</param>
        /// <param name="targetDB">The target database.</param>
        public static void CopyAll<T>(this IEnumerable<T> entities, IDatabaseInfo targetDB) {
            if (typeof(ISessionInfo).IsAssignableFrom(typeof(T))) {
                foreach (ISessionInfo session in entities) {
                    session.Database.Controller.CopySession(session, targetDB);
                }
            } else if (typeof(IGridInfo).IsAssignableFrom(typeof(T))) {
                foreach (IGridInfo grid in entities) {
                    grid.Database.Controller.CopyGrid(grid, targetDB);
                }
            } else {
                throw new ArgumentException(@"Wrong type of entities in this sequence.
                    As of now, you can only copy sessions or grids.");
            }
        }

        /// <summary>
        /// Moves all entities inside a collection to another database,
        /// provided that they are database entities that can be moves,
        /// i.e. ISessionInfo as of now.
        /// </summary>
        /// <typeparam name="T">The type of the entities.</typeparam>
        /// <param name="entities">The entities to be copied.</param>
        /// <param name="targetDB">The target database.</param>
        public static void MoveAll<T>(this IEnumerable<T> entities, IDatabaseInfo targetDB) {
            if (typeof(ISessionInfo).IsAssignableFrom(typeof(T))) {
                foreach (ISessionInfo session in entities) {
                    session.Database.Controller.MoveSession(session, targetDB);
                }
            } else {
                throw new ArgumentException(@"Wrong type of entities in this sequence.
                    As of now, you can only move sessions.");
            }
        }

        /// <summary>
        /// Returns the second entity in a collection.
        /// </summary>
        /// <typeparam name="T">The type of the entities.</typeparam>
        /// <param name="entities">The sequence of entities.</param>
        /// <returns>The second entity.</returns>
        public static T Second<T>(this IEnumerable<T> entities) {
            return entities.Pick(1);
        }

        /// <summary>
        /// Returns the third entity in a collection.
        /// </summary>
        /// <typeparam name="T">The type of the entities.</typeparam>
        /// <param name="entities">The sequence of entities.</param>
        /// <returns>The third entity.</returns>
        public static T Third<T>(this IEnumerable<T> entities) {
            return entities.Pick(2);
        }

        /// <summary>
        /// Returns the n-th entity in a collection.
        /// </summary>
        /// <typeparam name="T">The type of the entities.</typeparam>
        /// <param name="entities">The sequence of entities.</param>
        /// <param name="index">The (zero-based) index of the entity.</param>
        /// <returns>The n-th entity.</returns>
        public static T Pick<T>(this IEnumerable<T> entities, int index) {
            return entities.ElementAt(index);
        }

        /// <summary>
        /// Returns the elements at the given <paramref name="indices"/> from
        /// the given <paramref name="entities"/>.
        /// </summary>
        /// <typeparam name="T">The type of the entities.</typeparam>
        /// <param name="entities">The sequence of entities.</param>
        /// <param name="indices">
        /// The (zero-based) sequence of indices of the entity.
        /// </param>
        /// <returns>The elements at the given <paramref name="indices"/></returns>
        public static IEnumerable<T> Pick<T>(this IEnumerable<T> entities, params int[] indices) {
            return indices.Select(i => entities.ElementAt(i));
        }

        /// <summary>
        /// Variant of <see cref="Enumerable.Take"/> with a stride.
        /// </summary>
        /// <param name="timesteps"></param>
        /// <param name="noOfTimesteps"></param>
        /// <param name="stride"></param>
        /// <returns>
        /// Every <paramref name="stride"/>-th of the next
        /// <paramref name="noOfTimesteps"/> steps in
        /// <paramref name="timesteps"/>.
        /// </returns>
        public static IEnumerable<ITimestepInfo> Take(this IEnumerable<ITimestepInfo> timesteps, int noOfTimesteps, int stride) {
            return timesteps.Take(noOfTimesteps).Where((t, i) => i % stride == 0);
        }
    }
}
