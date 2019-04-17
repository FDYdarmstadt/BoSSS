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
using System.Diagnostics;

namespace BoSSS.Platform {

    /// <summary>
    /// A custom variant of the Lazy{T} class where the stored value is 
    /// automatically re-evaluated if it is no longer up-to-date.
    /// </summary>
    /// <typeparam name="T">
    /// The type of the lazy object
    /// </typeparam>
    public class ExpirableLazy<T> where T : class {

        /// <summary>
        /// Value of the object, if it has been created yet
        /// </summary>
        private T m_Value = null;

        /// <summary>
        /// The function for the creation of the value upon access.
        /// </summary>
        private Func<T> valueFun;

        /// <summary>
        /// The function that determines whether the represented value is
        /// still up to date.
        /// </summary>
        private Func<T, bool> isUpToDateFun;

        /// <summary>
        /// Constructs a wrapper for a lazy object.
        /// </summary>
        /// <param name="valueFun">
        /// The function for the creation of the value upon access.
        /// </param>
        /// <param name="isUpToDateFun">
        /// The function that determines whether the represented value is
        /// still up to date. If not, <paramref name="valueFun"/> will be
        /// called again in order to update the represented
        /// <see cref="Value"/>
        /// </param>
        public ExpirableLazy(Func<T> valueFun, Func<T, bool> isUpToDateFun) {
            this.valueFun = valueFun;
            this.isUpToDateFun = isUpToDateFun;
        }

        /// <summary>
        /// Returns true if the value has already been created.
        /// </summary>
        public bool IsValueCreated {
            get {
                return m_Value != null;
            }
        }

        /// <summary>
        /// Cached access to the lazy object.
        /// </summary>
        public T Value {
            get {
                if (m_Value == null) {


                    m_Value = valueFun();

                } else {
                    bool UpToDate = isUpToDateFun(m_Value);


                    if (!UpToDate) {
                        m_Value = valueFun();
                    }
                }

                return m_Value;
            }
        }
    }
}
