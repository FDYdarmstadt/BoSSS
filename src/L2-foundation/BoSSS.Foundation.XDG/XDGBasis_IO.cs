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
using BoSSS.Foundation.IO;

namespace BoSSS.Foundation.XDG {

    partial class XDGBasis {

        /// <summary>
        /// This class contains all necessary information to recreate an
        /// <see cref="XDGBasis"/>
        /// </summary>
        [Serializable]
        [DataContract]
        public class XDGBasisInitializer : BasisInitializer {

            /// <summary>
            /// Initializer of the tracker the represented basis depends on.
            /// </summary>
            [DataMember]
            public LevelSetTracker.LevelSetTrackerInitializer TrackerInitializer;

            /// <summary>
            /// See <see cref="Basis.BasisInitializer"/>
            /// </summary>
            /// <param name="c"></param>
            /// <returns></returns>
            public override Basis Initialize(IInitializationContext c) {
                XDGBasis xdgbasis;
                if (c.TryGetValue(this, out xdgbasis))
                    return xdgbasis;

                if (c.GridData == null)
                    throw new ArgumentException();
                if (!c.GridData.GridID.Equals(base.GridGuid))
                    throw new ArgumentException("Wrong grid.");

                var lsTrk = TrackerInitializer.Initialize(c);
                XDGBasis xb = new XDGBasis(lsTrk, base.Degree);
                myInstance = xb;
                c.Add(this, xb);
                return xb;
            }

            [NonSerialized]
            private XDGBasis myInstance;

            /// <summary>
            /// Compares the given object <paramref name="other"/> with respect
            /// to the 
            /// <see cref="BasisInitializer.GridGuid"/> and the 
            /// <see cref="BasisInitializer.Degree"/> and the 
            /// <see cref="TrackerInitializer"/>.
            /// </summary>
            public override bool Equals(Initializer<Basis> other) {
                XDGBasisInitializer initializer = other as XDGBasisInitializer;
                if(initializer == null)
                    return false;
                if(!initializer.GridGuid.Equals(this.GridGuid))
                    return false;
                if(!initializer.TrackerInitializer.Equals(this.TrackerInitializer))
                    return false;
                if(initializer.Degree != this.Degree)
                    return false;

                return true;
            }

            /// <summary>
            /// Computes a hash code based on 
            /// <see cref="GridGuid"/> and
            /// <see cref="Degree"/> and 
            /// <see cref="TrackerInitializer"/>
            /// </summary>
            public override int GetHashCode() {
                // http://stackoverflow.com/questions/1646807/quick-and-simple-hash-code-combinations
                int hash =  19937; // a prime number
                hash += 999331 * GridGuid.GetHashCode();
                hash += 999331 * base.Degree;
                hash += 999331 * this.TrackerInitializer.GetHashCode();
                
                return hash;
            }
        }

        /// <summary>
        /// Initializer of this basis object
        /// </summary>
        [DataMember]
        private XDGBasisInitializer m_Initializer;

        /// <summary>
        /// To support IO-architecture, NOT for direct user interaction. Note
        /// that it is essential that this member always returns the SAME
        /// object (reference-equals)!
        /// </summary>
        public override BasisInitializer Initializer {
            get {
                if (m_Initializer == null) {
                    m_Initializer = new XDGBasisInitializer() {
                        Degree = this.Degree,
                        GridGuid = this.GridDat.GridID,
                        TrackerInitializer = this.Tracker.Initializer
                    };
                }
                return m_Initializer;
            }
        }
    }
}
