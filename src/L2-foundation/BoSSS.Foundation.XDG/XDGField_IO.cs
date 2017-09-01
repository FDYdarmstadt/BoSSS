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
using System.Runtime.Serialization;
using System.Text;
using BoSSS.Foundation.IO;
using ilPSP;

namespace BoSSS.Foundation.XDG {

    public partial class XDGField {

        /// <summary>
        /// Specialized initializer for <see cref="XDGField"/>s
        /// </summary>
        [Serializable]
        public class XDGFieldInitializer : FieldInitializer {

            /// <summary>
            /// <see cref="DGField.FieldInitializer"/>
            /// </summary>
            /// <param name="c"></param>
            /// <returns></returns>
            public override DGField Initialize(IInitializationContext c) {
                DGField sff;
                if (c.TryGetValue(this, out sff))
                    return sff;

                var Basis = (XDGBasis)(base.BasisInfo.Initialize(c));
                XDGField f = new XDGField(Basis, this.Identification);
                myInstance = f;
                c.Add(this, f);
                return f;
            }

            [NonSerialized]
            private XDGField myInstance;

            /// <summary>
            /// Compares the given object <paramref name="other"/> with respect
            /// to the 
            /// <see cref="XDGFieldInitializer.Identification"/> and the 
            /// <see cref="FieldInitializer.BasisInfo"/>.
            /// </summary>
            /// <returns></returns>
            public override bool Equals(Initializer<DGField> other) {
                XDGFieldInitializer initializer = other as XDGFieldInitializer;
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
            /// <see cref="FieldInitializer.Identification"/> and 
            /// <see cref="FieldInitializer.BasisInfo"/>.
            /// </summary>
            public override int GetHashCode() {
                // http://stackoverflow.com/questions/1646807/quick-and-simple-hash-code-combinations
                int hash =  23; // a prime number
                hash += 113 * this.Identification.GetHashCode();
                hash += 113 * this.BasisInfo.GetHashCode();

                return hash;
            }
        }

        XDGFieldInitializer m_Initializer;

        /// <summary>
        /// To support IO-architecture, NOT for direct user interaction. Note
        /// that it is essential that this member always returns the SAME
        /// object (reference-equals)!
        /// </summary>
        public override DGField.FieldInitializer Initializer {
            get {
                if (m_Initializer == null) {
                    m_Initializer = new XDGFieldInitializer() {
                        BasisInfo = this.Basis.Initializer,
                        Identification = this.Identification
                    };
                }
                return m_Initializer;
            }
        }

        /// <summary>
        /// <see cref="DGField.FieldInitializer"/>
        /// </summary>
        /// <param name="tsi"></param>
        /// <param name="data"></param>
        /// <param name="loadedObjects"></param>
        public override void LoadData(ITimestepInfo tsi, IList<CellFieldDataSet> data, HashSet<object> loadedObjects) {
            if (loadedObjects.Contains(this))
                return;

            this.Basis.Tracker.LoadData(tsi, data, loadedObjects);

            int MyIndex = tsi.FieldInitializers.IndexOf(
                this.Initializer, (a, b) => a.Identification.Equals(b.Identification));
            XDGFieldInitializer myInfo = (XDGFieldInitializer)tsi.FieldInitializers.Single(info =>
                info.Identification.Equals(this.Identification));

            if (this.Basis.Degree == myInfo.BasisInfo.Degree) {
                XDGField dis = this;
                LoadCoordinates(data, MyIndex, dis);
            } else {
                XDGField Temp = new XDGField(new XDGBasis(this.Basis.Tracker, myInfo.BasisInfo.Degree));
                LoadCoordinates(data, MyIndex, Temp);
                this.Clear();
                this.AccLaidBack(1.0, Temp);
            }
            loadedObjects.Add(this);
        }

        /// <summary>
        /// Loads the XDG coordinates for the given <paramref name="field"/>
        /// from the given block of <paramref name="data"/>.
        /// </summary>
        /// <param name="data"></param>
        /// <param name="MyIndex"></param>
        /// <param name="field"></param>
        private static void LoadCoordinates(IList<CellFieldDataSet> data, int MyIndex, XDGField field) {
            int J = field.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            for (int j = 0; j < J; j++) {
                //var coords_j = data[j].DGCoordinateData[MyIndex];
                var coords_j = data[j].GetDGCoordinates(MyIndex);

                //if (coords_j.Length != field.Basis.GetLength(j))
                if (coords_j.Length != field.Basis.GetLength(j)) {
                    throw new Exception();
                    //Console.WriteLine("Bullshit in cell {0}", j);
                    //field.Coordinates.ClearRow(j);
                } else {
                    for (int n = 0; n < coords_j.Length; n++)
                        field.Coordinates[j, n] = coords_j[n];
                }
            }
        }

        /// <summary>
        /// Depends on the level set(s)
        /// </summary>
        /// <returns></returns>
        public override IEnumerable<DGField> ReportDependentFields() {
            return this.Basis.Tracker.LevelSets.Select(ls => (LevelSet)ls);
        }
    }
}
