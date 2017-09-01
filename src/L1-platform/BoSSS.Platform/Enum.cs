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

namespace BoSSS.Platform {

    /// <summary>
    /// Utility methods for <see cref="Enum"/>s exploiting generics.
    /// </summary>
    /// <typeparam name="EnumType">
    /// The type of the considered enum
    /// </typeparam>
    public static class Enum<EnumType> where EnumType : struct {

        /// <summary>
        /// Parses the given string value into an enum of type
        /// <typeparamref name="EnumType"/>. The case of the given value
        /// <paramref name="value"/> is considered relevant.
        /// </summary>
        /// <param name="value">The value to be parsed</param>
        /// <returns>The parsed enum</returns>
        public static EnumType Parse(string value) {
            return Parse(value, false);
        }

        /// <summary>
        /// Parses the given string value into an enum of type
        /// <typeparamref name="EnumType"/>. Depending on
        /// <paramref name="ignoreCase"/>, the case of
        /// <paramref name="value"/> is or is not considered relevant.
        /// </summary>
        /// <param name="value">The value to be parsed</param>
        /// <param name="ignoreCase">
        /// <see cref="Enum.Parse(Type, string, bool) "/>
        /// </param>
        /// <returns>The parsed enum</returns>
        public static EnumType Parse(string value, bool ignoreCase) {
            CheckIsEnum();
            return (EnumType)Enum.Parse(typeof(EnumType), value, ignoreCase);
        }

        /// <summary>
        /// Tries to parse the given string value (case insensitive comparison)
        /// and returns the default value of <typeparamref name="EnumType"/> in
        /// case this was not successful.
        /// </summary>
        /// <param name="value">The value to be parsed</param>
        /// <returns>
        /// The enum value corresponding to <paramref name="value"/> if it
        /// exists. Otherwise, default(<typeparamref name="EnumType"/>) is
        /// returned.
        /// </returns>
        public static EnumType ParseOrDefault(string value) {
            EnumType result;
            TryParse(value, out result);
            return result;
        }

        /// <summary>
        /// Tries to parse the given string value (case insensitive comparison)
        /// </summary>
        /// <param name="value">The value to be parsed</param>
        /// <param name="result">
        /// On exit: the enum value corresponding to <paramref name="value"/>
        /// if it exists. Otherwise, default(<typeparamref name="EnumType"/>).
        /// </param>
        /// <returns>
        /// True, if the parsing of <paramref name="value"/> was successful;
        /// false otherwise.
        /// </returns>
        public static bool TryParse(string value, out EnumType result) {
            foreach (EnumType enumValue in GetValues()) {
                if (enumValue.ToString().Equals(value, StringComparison.InvariantCultureIgnoreCase)) {
                    result = enumValue;
                    return true;
                }
            }

            result = default(EnumType);
            return false;
        }

        /// <summary>
        /// Type-safe variant of <see cref="Enum.GetValues"/>.
        /// </summary>
        /// <returns>
        /// The list of constants defined in <typeparamref name="EnumType"/>.
        /// </returns>
        public static IEnumerable<EnumType> GetValues() {
            IList<EnumType> list = new List<EnumType>();
            foreach (object value in Enum.GetValues(typeof(EnumType))) {
                list.Add((EnumType)value);
            }
            return list;
        }

        /// <summary>
        /// Checks whether the passed type parameter
        /// <typeparamref name="EnumType"/> really is an enum since this
        /// generic type constraint is not supported by C#.
        /// </summary>
        private static void CheckIsEnum() {
            if (!typeof(EnumType).IsEnum) {
                throw new Exception("Type parameter " + typeof(EnumType).Name + " must be an enum");
            }
        }
    }
}
