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
using BoSSS.Foundation.IO;
using BoSSS.Platform;

namespace BoSSS.Foundation.XDG {

    public partial class LevelSetTracker {

        /// <summary>
        /// Initializer for a <see cref="LevelSetTracker"/> during IO
        /// operations
        /// </summary>
        [Serializable]
        [DataContract]
        public class LevelSetTrackerInitializer : Initializer<LevelSetTracker> {

            /// <summary>
            /// An instance of the represented tracker (if already constructed
            /// via <see cref="Initialize"/>).
            /// </summary>
            [NonSerialized]
            private LevelSetTracker instance;

            /// <summary>
            /// See <see cref="LevelSetTracker.LevelSets"/>
            /// </summary>
            [DataMember]
            public LevelSet.LevelSetInitializer[] LevelSets;

            /// <summary>
            /// See <see cref="LevelSetTracker.NearRegionWidth"/>
            /// </summary>
            [DataMember]
            internal int NearRegionWidth;

            /// <summary>
            /// See <see cref="LevelSetTracker.CutCellQuadratureType"/>
            /// </summary>
            [DataMember]
            public XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;

            /// <summary>
            /// See <see cref="LevelSetTracker.SpeciesTable"/>
            /// </summary>
            [DataMember]
            internal Array SpeciesTable;

            /// <summary>
            /// Initializes the level set tracker.
            /// </summary>
            /// <param name="c"></param>
            /// <returns></returns>
            public override LevelSetTracker Initialize(IInitializationContext c) {
                LevelSetTracker lstrk;
                if (c.TryGetValue(this, out lstrk))
                    return lstrk;

                LevelSet[] LS = new LevelSet[this.LevelSets.Length];
                for (int i = 0; i < LS.Length; i++) {
                    LS[i] = (LevelSet)this.LevelSets[i].Initialize(c);
                }

                var lsTrk = new LevelSetTracker(c.GridData, this.CutCellQuadratureType, this.NearRegionWidth, this.SpeciesTable, LS);
                instance = lsTrk;
                c.Add(this, lsTrk);
                return lsTrk;
            }

            /// <summary>
            /// Compares the given object <paramref name="other"/> with this one.
            /// </summary>
            public override bool Equals(Initializer<LevelSetTracker> other) {
                LevelSetTrackerInitializer initializer = other as LevelSetTrackerInitializer;
                if (initializer == null)
                    return false;
                if (this.LevelSets.Length != initializer.LevelSets.Length)
                    return false;
                for (int i = 0; i < this.LevelSets.Length; i++) {
                    if (!this.LevelSets[i].Equals(initializer.LevelSets[i]))
                        return false;
                }
                if(this.CutCellQuadratureType != initializer.CutCellQuadratureType)
                    return false;
                if (this.NearRegionWidth != initializer.NearRegionWidth)
                    return false;
                if (this.SpeciesTable.Length != initializer.SpeciesTable.Length)
                    return false;
                if (this.SpeciesTable.Rank != initializer.SpeciesTable.Rank)
                    return false;
                int R = this.SpeciesTable.Rank;
                int L = this.SpeciesTable.Length;
                int[] LT = new int[R];
                for (int r = 0; r < R; r++) {
                    if (this.SpeciesTable.GetLength(r) != initializer.SpeciesTable.GetLength(r))
                        return false;
                }

                LT[R - 1] = 1;
                for (int r = R - 2; r >= 0; r--)
                    LT[r] = this.SpeciesTable.GetLength(r + 1) * LT[r + 1];

                int[] idx = new int[R];
                for (int l = 0; l < L; l++) {
                    int ll = l;
                    for (int r = 0; r < R; r++) {
                        idx[r] = ll / LT[r];
                        ll -= idx[r] * LT[r];
                    }

                    string sa = (string)this.SpeciesTable.GetValue(idx);
                    string sb = (string)initializer.SpeciesTable.GetValue(idx);

                    if (sa != sb)
                        return false;
                }

                return true;
            }

            /// <summary>
            /// Computes a hash code based on 
            /// <see cref="LevelSets"/> and
            /// <see cref="NearRegionWidth"/>.
            /// </summary>
            public override int GetHashCode() {
                // http://stackoverflow.com/questions/1646807/quick-and-simple-hash-code-combinations
                int hash = 311; // a prime number
                foreach (var ls in this.LevelSets)
                    hash += 1931 * ls.GetHashCode();
                hash += 1931 * this.NearRegionWidth;

                return hash;
            }
        }

        /// <summary>
        /// Initializer of this object, if any.
        /// </summary>
        private LevelSetTrackerInitializer m_Initializer;

        /// <summary>
        /// To support IO-architecture, NOT for direct user interaction. Note
        /// that it is essential that this member always returns the SAME
        /// object (reference-equals)!
        /// </summary>
        public LevelSetTrackerInitializer Initializer {
            get {
                if (m_Initializer == null) {
                    m_Initializer = new LevelSetTrackerInitializer() {
                        LevelSets = this.LevelSets.Select(ls =>
                            (LevelSet.LevelSetInitializer)ls.As<LevelSet>().Initializer).ToArray(),
                        NearRegionWidth = this.NearRegionWidth,
                        SpeciesTable = this.m_SpeciesTable,
                        CutCellQuadratureType = this.CutCellQuadratureType
                    };
                }
                return m_Initializer;
            }
        }

        /// <summary>
        /// Makes each level set load the appropriate <paramref name="data"/>
        /// for the given <paramref name="tsi"/>.
        /// </summary>
        /// <param name="tsi">
        /// Information about the time-step
        /// </param>
        /// <param name="data">
        /// Data to be loaded
        /// </param>
        /// <param name="loadedObjects">
        /// Cache for already loaded objects
        /// </param>
        /// <remarks>
        /// Causes an update of the level set tracker!
        /// </remarks>
        public void LoadData(ITimestepInfo tsi, IList<CellFieldDataSet> data, HashSet<object> loadedObjects) {
            if (loadedObjects.Contains(this))
                return;

            foreach (var Ls in this.LevelSets) {
                Ls.As<LevelSet>().LoadData(tsi, data, loadedObjects);
            }
            this.UpdateTracker();
            loadedObjects.Add(this);
        }
    }
}
