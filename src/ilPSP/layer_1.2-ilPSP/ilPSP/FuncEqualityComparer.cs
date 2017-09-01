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

namespace ilPSP {

    

    /// <summary>
    /// An equality comparer that can be build from a lambda-expression.
    /// </summary>
    public class FuncEqualityComparer<T> : IEqualityComparer<T> {

        /// <summary>
        /// The comparer to be used
        /// </summary>
        private readonly Func<T, T, bool> _comparer;

        /// <summary>
        /// An optional hashing function.
        /// </summary>
        private readonly Func<T, int> _hash;

        /// <summary>
        /// Constructs an equality comparer without any hash function
        /// </summary>
        /// <param name="comparer">An equality comparison function</param>
        public FuncEqualityComparer(Func<T, T, bool> comparer)
            : this(comparer, t => 0)
        {
        }

        /// <summary>
        /// Constructs an equality comparer with a selected hash function
        /// </summary>
        /// <param name="comparer">An equality comparison function</param>
        /// <param name="hash">A hash function</param>
        public FuncEqualityComparer(Func<T, T, bool> comparer, Func<T, int> hash) {
            _comparer = comparer;
            _hash = hash;
        }

        /// <summary>
        /// Equality
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public bool Equals(T x, T y) {
            return _comparer(x, y);
        }

        /// <summary>
        /// Hash code
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public int GetHashCode(T obj) {
            return _hash(obj);
        }
    }

    /// <summary>
    /// some extension functions
    /// </summary>
    public static class FuncEqualityComparerExtensions {

        /// <summary>
        /// Alternative to the constructor of <see cref="FuncEqualityComparer{T}"/>
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="equality"></param>
        /// <returns></returns>
        public static FuncEqualityComparer<T> ToEqualityComparer<T>(this Func<T, T, bool> equality) {
            return new FuncEqualityComparer<T>(equality);
        }
    }
}
