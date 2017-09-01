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
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Text;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Extension methods for <see cref="Type"/>.
    /// </summary>
    public static class TypeExtensions {
        
        /// <summary>
        /// Retrieves the extension methods for <paramref name="type"/> in the
        /// given <paramref name="assembly"/>.
        /// </summary>
        /// <param name="type">
        /// The type of interest.
        /// </param>
        /// <param name="assembly">
        /// The assembly potentially containing extension methods.
        /// </param>
        /// <returns>
        /// A (potentially incomplete) list of extension methods for the given
        /// <paramref name="type"/>.
        /// </returns>
        /// <remarks>
        /// Unfortunately, finding the full list of extension methods via
        /// reflection is not that easy since it is hard to determine whether
        /// a given type satisfies the generic restriction of a method (i.e.,
        /// <see cref="Type.IsAssignableFrom"/> does not cover this case). At
        /// the moment, the corresponding extension methods will not be found.
        /// </remarks>
        public static IEnumerable<MethodInfo> GetExtensionMethods(this Type type, Assembly assembly) {
            BindingFlags flags = BindingFlags.Static | BindingFlags.Public | BindingFlags.NonPublic;
            return assembly.GetTypes().SelectMany(t =>
                t.GetMethods(flags).Where(m =>
                    m.IsDefined(typeof(ExtensionAttribute), false)
                    && m.GetParameters().First().ParameterType.IsAssignableFrom(type)));
        }

        /// <summary>
        /// For a generic type <paramref name="t"/> that represents
        /// someType{T}, retrieves the name ('someType')
        /// </summary>
        /// <param name="t">The considered type</param>
        /// <param name="fullyQualified">
        /// If true, the type names will be prefixed with its full namespace.
        /// </param>
        /// <returns>
        /// The name of the generic type.
        /// </returns>
        public static string GetBaseName(this Type t, bool fullyQualified = false) {
            if (!t.IsGenericType) {
                return t.Name;
            }

            Type type = t.GetGenericTypeDefinition();
            string genericName = fullyQualified ? t.FullName : t.Name;
            return genericName.Substring(0, genericName.IndexOf('`'));
        }

        /// <summary>
        /// Recursively assembles the full, human-readable type name of a
        /// (generic) class. For example, returns 'IEnumerable{SomeClass}'
        /// or 'KeyValuePair{SomeClass, SomeOtherClass}' instead of the cryptic
        /// names used by the CLR.
        /// </summary>
        /// <param name="t">The type whose name is of interest</param>
        /// <param name="fullyQualified">
        /// If true, all appearing type names will be prefixed with their full
        /// namespace.
        /// </param>
        /// <returns>
        /// The human-readable name of <paramref name="t"/>.
        /// </returns>
        public static string GetGenericTypeName(this Type t, bool fullyQualified = false) {
            if (!t.IsGenericType) {
                return fullyQualified ? t.FullName : t.Name;
            }

            StringBuilder name = new StringBuilder(t.GetBaseName(fullyQualified));
            var typeArgs = t.GetGenericArguments();
            name.Append("{" + typeArgs.First().GetGenericTypeName(fullyQualified));
            foreach (var typeArg in typeArgs.Skip(1)) {
                name.Append("," + typeArg.GetGenericTypeName(fullyQualified));
            }
            name.Append("}");

            return name.ToString();
        }
    }
}
