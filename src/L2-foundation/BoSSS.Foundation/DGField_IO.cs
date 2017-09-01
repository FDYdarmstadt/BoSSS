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
using System.Diagnostics;
using BoSSS.Foundation.IO;
using ilPSP;

namespace BoSSS.Foundation {

    public partial class DGField {

        /// <summary>
        /// To support IO-architecture, NOT for direct user interaction. 
        /// Note that it is essential that this member always returns the SAME
        /// object (reference-equals)!
        /// </summary>
        public abstract FieldInitializer Initializer {
            get;
        }

        /// <summary>
        /// Represents a set of information about a field that sufficient for
        /// its recreation, i.e. during deserialization.
        /// </summary>
        [Serializable]
        public abstract class FieldInitializer : Initializer<DGField> {

            /// <summary>
            /// Information about the basis.
            /// </summary>
            public Basis.BasisInitializer BasisInfo;

            /// <summary>
            /// The name of the field
            /// </summary>
            public string Identification;

            /// <summary>
            /// Must be overridden.
            /// </summary>
            public override bool Equals(Initializer<DGField> other) {
                throw new NotImplementedException("method must be overridden.");
            }

            /// <summary>
            /// Must be overridden.
            /// </summary>
            public override int GetHashCode() {
                throw new NotImplementedException("method must be overridden.");
            }
        }

        /// <summary>
        /// Supports the loading of DG fields from database, not intended for
        /// direct user interaction.
        /// </summary>
        public virtual void LoadData(ITimestepInfo tsi, IList<CellFieldDataSet> data, HashSet<object> loadedObjects) {
            if (loadedObjects.Contains(this))
                return;

            if (this.Identification == null || this.Identification.Length <= 0)
                throw new NotSupportedException("unable to load a timestep into unnamed fields.");

            int MyIndex = tsi.FieldInitializers.IndexOf(
                this.Initializer, (a, b) => a.Identification.Equals(b.Identification));

            if (MyIndex >= 0) {
                int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                for (int j = 0; j < J; j++) {
                    //double[] coords_j = data[j].DGCoordinateData[MyIndex];
                    //double[] coords_j = data[j].DGCoordinateData[MyIndex].Data;
                    double[] coords_j = data[j].GetDGCoordinates(MyIndex);

                    Debug.Assert(data[j].GlobalID == this.GridDat.iLogicalCells.GetGlobalID(j));

                    int Nt = this.Basis.GetLength(j);
                    int N = Math.Min(coords_j.Length, Nt);

                    int n = 0;
                    for (; n < N; n++) {
                        this.Coordinates[j, n] = coords_j[n];
                    }
                    for (; n < Nt; n++) {
                        this.Coordinates[j, n] = 0.0;
                    }
                }
            } else {
                Console.WriteLine("Unable to load field '{0}'; initializing with zeros.", this.Identification);
            }
            loadedObjects.Add(this);
        }

        /// <summary>
        /// To support the database loading architecture, not for direct user
        /// interaction. Example use case: an XDG field reports the level sets
        /// it depends on.
        /// </summary>
        public virtual IEnumerable<DGField> ReportDependentFields() {
            return new DGField[0];
        }
    }
}
