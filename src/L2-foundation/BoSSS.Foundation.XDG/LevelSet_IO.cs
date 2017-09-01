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
using BoSSS.Foundation.IO;

namespace BoSSS.Foundation.XDG {

    public partial class LevelSet {

        /// <summary>
        /// See <see cref="Initializer"/>
        /// </summary>
        LevelSetInitializer m_Initializer;

        /// <summary>
        /// To support IO-architecture, NOT for direct user interaction. Note
        /// that it is essential that this member always returns the SAME
        /// object (reference-equals)!
        /// </summary>
        public override DGField.FieldInitializer Initializer {
            get {
                if (m_Initializer == null) {
                    m_Initializer = new LevelSetInitializer() {
                        BasisInfo = this.Basis.Initializer,
                        Identification = this.Identification
                    };
                }
                return m_Initializer;
            }
        }

        /// <summary>
        /// Specialized initializer for level set fields.
        /// </summary>
        [Serializable]
        public class LevelSetInitializer : SinglePhaseFieldInitializer {

            /// <summary>
            /// <see cref="IInitializer{T}.Initialize"/>
            /// </summary>
            /// <param name="c"></param>
            /// <returns></returns>
            public override DGField Initialize(IInitializationContext c) {
                DGField sff;
                if (c.TryGetValue(this, out sff))
                    return sff;

                var Basis = base.BasisInfo.Initialize(c);
                LevelSet ls = new LevelSet(Basis, this.Identification);
                base.myInstance = ls;
                c.Add(this, ls);
                return ls;
            }

            /// <summary>
            /// Compares the given object <paramref name="other"/> with respect
            /// to the 
            /// <see cref="DGField.FieldInitializer.Identification"/> and the 
            /// <see cref="DGField.FieldInitializer.BasisInfo"/>.
            /// </summary>
            /// <returns></returns>
            public override bool Equals(Initializer<DGField> other) {
                LevelSetInitializer initializer = other as LevelSetInitializer;
                if (initializer == null)
                    return false;
                if (!base.BasisInfo.Equals(initializer.BasisInfo))
                    return false;
                if (!base.Identification.Equals(initializer.Identification))
                    return false;

                return true;
            }

            /// <summary>
            /// Computes a hash code based on 
            /// <see cref="DGField.FieldInitializer.Identification"/> and 
            /// <see cref="DGField.FieldInitializer.BasisInfo"/>.
            /// </summary>
            public override int GetHashCode() {
                // http://stackoverflow.com/questions/1646807/quick-and-simple-hash-code-combinations
                int hash = 373; // a prime number
                hash += 719 * this.Identification.GetHashCode();
                hash += 719 * this.BasisInfo.GetHashCode();

                return hash;
            }
        }
    }
}
