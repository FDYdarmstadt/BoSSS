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

namespace ilPSP {
    
    
    /// <summary>
    /// creates a comparer-object from a delegate
    /// </summary>
    public class FuncComparer<T> : IComparer<T> {

        Func<T, T, int> m_comp;

        /// <summary>
        /// Constructs a comparer
        /// </summary>
        /// <param name="comparer">An comparison function</param>
        public FuncComparer(Func<T, T, int> comparer)
        {
            this.m_comp = comparer;
        }

        /// <summary>
        /// comparison implementation
        /// </summary>
        public int Compare(T x, T y) {
            return m_comp(x,y);
        }
    }


    /// <summary>
    /// some extension functions
    /// </summary>
    public static class FuncComparerExtensions {

        /// <summary>
        /// Alternative to the constructor of <see cref="FuncComparer{T}"/>
        /// </summary>
        public static FuncComparer<T> ToComparer<T>(this Func<T, T, int> c) {
            return new FuncComparer<T>(c);
        }
    }
}
