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
            /// Compares the given object <paramref name="other"/> 
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
            /// Computes a hash code
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
        /// stores XDG coordinates species-wise and retains backward compatibility using a mysterious magic header
        /// </summary>
        public override double[] SerializeDGcoords(int j) {

            int Ndg = this.Basis.NonX_Basis.GetLength(j);
            var trk = this.Basis.Tracker;
            int NoOfSpc = trk.Regions.GetNoOfSpecies(j);
            


            double[] Ret = new double[4 + NoOfSpc*(Ndg +1)];

            // write the magic header
            Ret[0] = double.NegativeInfinity;
            Ret[1] = double.MaxValue;
            Ret[2] = double.NaN;
            Ret[3] = Ndg;

            int Ptr = 4;
            for(int iSpc = 0; iSpc < NoOfSpc; iSpc++) {
                int n0 = iSpc * Ndg;
                SpeciesId spc = trk.Regions.GetSpeciesIdFromIndex(j, iSpc);
                if(trk.Regions.IsSpeciesPresentInCell(spc, j)) {
                    Ret[Ptr] = spc.cntnt; Ptr++;
                    for(int n = 0; n < Ndg; n++) {
                        Ret[Ptr] = this.Coordinates[j, n + n0];
                        Ptr++;
                    }
                }
            }

            if(Ptr < Ret.Length)
                Array.Resize(ref Ret, Ptr);

            return Ret;
        }

        /// <summary>
        /// loads XDG coordinates species-wise and retains backward compatibility using a mysterious magic header
        /// </summary>
        public override void DeserializeDGcoords(int j, double[] coords_j) {
            if(coords_j.Length >= 3 
                && double.IsNegativeInfinity(coords_j[0]) && coords_j[1] == double.MaxValue && double.IsNaN(coords_j[2])) {
                // ++++++++++++++++++
                // magic header found
                // ++++++++++++++++++

                var trk = this.Basis.Tracker;

                int NdgStore = (int) coords_j[3];
                int NdgAct = this.Basis.NonX_Basis.GetLength(j);
                int Ndg = Math.Min(NdgAct, NdgStore);
                int Ptr = 4;

                while(Ptr < coords_j.Length) {

                    SpeciesId spc;
                    spc.cntnt = (int)coords_j[Ptr]; Ptr++;
                    int iSpc = trk.Regions.GetSpeciesIndex(spc, j);

                    int n0 = iSpc * NdgAct;
                    for(int n = 0; n < Ndg; n++) {
                        this.Coordinates[j, n0 + n] = coords_j[Ptr];
                        Ptr++;
                    }

                    for(int n = Ndg; n < NdgStore; n++)
                        Ptr++;
                    
                }

            } else {
                // +++++++++++++++++++++++++++
                // no magic header:
                // seems to be the old version
                // +++++++++++++++++++++++++++


                var trk = this.Basis.Tracker;

                int NoSpc = trk.Regions.GetNoOfSpecies(j);

                int NdgStore = coords_j.Length / NoSpc;
                int NdgAct = this.Basis.NonX_Basis.GetLength(j);
                int Nload = Math.Min(NdgAct, NdgStore);
                int Ptr = 0;

                for(int iSpc = 0; iSpc < NoSpc; iSpc++) {
                    int n0 = iSpc * NdgAct;
                    for(int n = 0; n < Nload; n++) {
                        this.Coordinates[j, n0 + n] = coords_j[Ptr];
                        Ptr++;
                    }

                    for(int n = Nload; n < NdgStore; n++)
                        Ptr++;
                }
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

            if (MyIndex >= 0) {
                int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                for (int j = 0; j < J; j++) {
                    //double[] coords_j = data[j].DGCoordinateData[MyIndex];
                    //double[] coords_j = data[j].DGCoordinateData[MyIndex].Data;
                    double[] coords_j = data[j].GetDGCoordinates(MyIndex);

                    Debug.Assert(data[j].GlobalID == this.GridDat.iLogicalCells.GetGlobalID(j));

                    this.DeserializeDGcoords(j, coords_j);
                }
            } else {
                Console.Error.WriteLine("Unable to load field '{0}'; initializing with zeros.", this.Identification);
                this.Coordinates.Clear();
            }


            /*
            if (this.Basis.Degree == myInfo.BasisInfo.Degree) {
                XDGField dis = this;
                LoadCoordinates(data, MyIndex, dis);
            } else {
                XDGField Temp = new XDGField(new XDGBasis(this.Basis.Tracker, myInfo.BasisInfo.Degree));
                LoadCoordinates(data, MyIndex, Temp);
                this.Clear();
                this.AccLaidBack(1.0, Temp);
            }*/

            loadedObjects.Add(this);
        }

        /*
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

                int cjL = coords_j.Length;
                //if (coords_j.Length != field.Basis.GetLength(j))
                if (cjL != field.Basis.GetLength(j)) {
                    //throw new Exception();
                    //Console.WriteLine("Bullshit in cell {0}", j);
                    //field.Coordinates.ClearRow(j);
                    ReducedRegionCode dummy;
                    int NoSpc = field.Basis.Tracker.Regions.GetNoOfSpecies(j, out dummy);
                    if (cjL == field.Basis.GetLength(j) / NoSpc) {
                        Console.WriteLine("Warning: field marked as cut-cell but not sufficient data to load");
                        for (int spc = 0; spc < NoSpc; spc++) {
                            for (int n = 0; n < cjL; n++)
                                field.Coordinates[j, (spc * cjL) + n] = coords_j[n];
                        }
                    } else {
                        throw new Exception();
                    }
                } else {
                    for (int n = 0; n < cjL; n++)
                        field.Coordinates[j, n] = coords_j[n];
                }
            }
        }
        */

        /// <summary>
        /// Depends on the level set(s)
        /// </summary>
        /// <returns></returns>
        public override IEnumerable<DGField> ReportDependentFields() {
            return this.Basis.Tracker.LevelSets.Select(ls => (LevelSet)ls);
        }
    }
}
