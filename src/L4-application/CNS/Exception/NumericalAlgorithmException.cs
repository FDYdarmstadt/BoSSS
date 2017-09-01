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
using System.Runtime.Serialization;

namespace CNS.Exception {

    /// <summary>
    /// Exception thrown if a numerical algorithm fails (e.g. because of
    /// unphsyical conditions)
    /// </summary>
    [Serializable]
    public class NumericalAlgorithmException : System.Exception {
        
        /// <summary>
        /// <see cref="Exception"/>
        /// </summary>
        public NumericalAlgorithmException()
            : base() {
        }

        /// <summary>
        /// <see cref="Exception"/>
        /// </summary>
        /// <param name="message"><see cref="Exception"/></param>
        public NumericalAlgorithmException(string message)
            : base(message) {
        }

        /// <summary>
        /// <see cref="Exception"/>
        /// </summary>
        /// <param name="message"><see cref="Exception"/></param>
        /// <param name="innerException"><see cref="Exception"/></param>
        public NumericalAlgorithmException(string message, System.Exception innerException)
            : base(message, innerException) {
        }

        /// <summary>
        /// <see cref="Exception"/>
        /// </summary>
        /// <param name="info"><see cref="Exception"/></param>
        /// <param name="context"><see cref="Exception"/></param>
        public NumericalAlgorithmException(SerializationInfo info, StreamingContext context)
            : base(info, context) {
        }
    }
}
