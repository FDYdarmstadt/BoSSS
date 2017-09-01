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
using ilPSP;

namespace BoSSS.Platform {

    /// <summary>
    /// Extension methods for <see cref="object"/>
    /// </summary>
    public static class ObjectExtensions {

        /// <summary>
        /// Provides a more convenient syntax for cast operations.
        /// </summary>
        /// <typeparam name="T">
        /// The type to be casted to
        /// </typeparam>
        /// <param name="obj">
        /// The object to be casted
        /// </param>
        /// <returns>The casted object</returns>
        public static T Cast<T>(this object obj) {
            return (T)obj;
        }

        /// <summary>
        /// Provides a more convenient syntax for 'as' operations.
        /// </summary>
        /// <typeparam name="T">
        /// The type to be converted to
        /// </typeparam>
        /// <param name="obj">
        /// The object to be converted to
        /// </param>
        /// <returns>The converted object</returns>
        public static T As<T>(this object obj) where T : class {
            return obj as T;
        }

        /// <summary>
        /// Returns the result of <paramref name="function"/>(<paramref name="obj"/>),
        /// if <paramref name="obj"/> is not null;
        /// otherwise, the default value of <typeparamref name="R"/>.
        /// </summary>
        public static R IfNotNull<T, R>(this T obj, Func<T, R> function)
            where T : class {
            if (obj == null) {
                return default(R);
            } else {
                return function(obj);
            }
        }

        /// <summary>
        /// Turns a single element into something that is enumerable (i.e., no
        /// array or something needs to be created by hand)
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="obj"></param>
        /// <returns>
        /// A sequence containing <paramref name="obj"/> only.
        /// </returns>
        public static IEnumerable<T> ToEnumerable<T>(this T obj) {
            yield return obj;
        } 
    }
}
