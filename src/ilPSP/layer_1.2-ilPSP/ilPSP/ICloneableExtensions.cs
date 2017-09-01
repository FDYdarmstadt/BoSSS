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

namespace ilPSP {

    /// <summary>
    /// Extension methods for the <see cref="ICloneable"/> interface.
    /// </summary>
    public static class ICloneableExtensions {

        /// <summary>
        /// Clones an object and casts it to the type of the original object
        /// at the same time.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the cloneable object and its clone
        /// </typeparam>
        /// <param name="obj">
        /// The object to be cloned
        /// </param>
        /// <returns>
        /// A clone of <paramref name="obj" /> (see
        /// <see cref="ICloneable.Clone" />)
        /// </returns>
        /// <exception cref="System.ArgumentException">
        /// Object must not be null
        /// </exception>
        public static T CloneAs<T>(this T obj) where T : ICloneable {
            if (obj == null) {
                throw new ArgumentException("Object must not be null", "obj");
            }

            return (T)obj.Clone();
        }

        /// <summary>
        /// Clones as or default.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the cloneable object and its clone
        /// </typeparam>
        /// <param name="obj">
        /// The object to be cloned
        /// </param>
        /// <returns>
        /// A clone of <paramref name="obj" /> (see
        /// <see cref="ICloneable.Clone" />) or the default value for the type
        /// <typeparamref name="T"/> if <paramref name="obj"/> is null.
        /// </returns>
        public static T CloneAsOrDefault<T>(this T obj) where T : ICloneable {
            if (obj == null) {
                return default(T);
            } else {
                return (T)obj.Clone();
            }
        }
    }
}
