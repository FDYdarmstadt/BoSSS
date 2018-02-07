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

namespace BoSSS.Solution.Control {

    /// <summary>
    /// Options for DG fields
    /// </summary>
    [Serializable]
    [DataContract]
    public class FieldOpts {

        /// <summary>
        /// DG polynomial degree; below zero denotes not specified.
        /// </summary>
        [DataMember]
        public int Degree = -1;

        /// <summary>
        /// see <see cref="SaveToDB"/>
        /// </summary>
        public enum SaveToDBOpt {
            /// <summary> 
            /// %
            /// </summary>
            TRUE,

            /// <summary>
            /// %
            /// </summary>
            FALSE,

            /// <summary>
            /// user has not specified anything in ctrl file
            /// </summary>
            unspecified
        }

        /// <summary>
        /// field state should be saved in control file
        /// </summary>
        [DataMember]
        public SaveToDBOpt SaveToDB = SaveToDBOpt.unspecified;

        /// <summary>
        /// %
        /// </summary>
        public override bool Equals(object obj) {
            FieldOpts f = obj as FieldOpts;
            if (f == null)
                return false;

            if (Degree != f.Degree)
                return false;
            if (SaveToDB != f.SaveToDB)
                return false;

            return true;
        }

        /// <summary>
        /// %
        /// </summary>
        public override int GetHashCode() {
            return Degree;
        }

    }
}
