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
using BoSSS.Foundation;

namespace BoSSS.Application.BoSSSpad {
    
    /// <summary>
    /// Extension methods for <see cref="DGField"/>
    /// </summary>
    public static class DGFieldExtensions {

        /// <summary>
        /// Retrieves the DG field with the given <paramref name="name"/>.
        /// </summary>
        /// <param name="fields">List to be searched</param>
        /// <param name="name">Name to be searched for</param>
        /// <returns>
        /// The single field within <paramref name="fields"/> with the given
        /// <paramref name="name"/>
        /// </returns>
        public static DGField Find(this IEnumerable<DGField> fields, string name) {
            return fields.Where(f => f.Identification == name).SingleOrDefault();
        }

        /// <summary>
        /// Retrieves the DG fields whose names start with the given
        /// <paramref name="name"/>.
        /// </summary>
        /// <param name="fields">List to be searched</param>
        /// <param name="name">Name to be searched for</param>
        /// <returns>
        /// The fields within <paramref name="fields"/> whose name starts with
        /// the given <paramref name="name"/>
        /// </returns>
        public static IEnumerable<DGField> FindAll(this IEnumerable<DGField> fields, string name) {
            return fields.Where(f => f.Identification.StartsWith(name));
        }
    }
}
