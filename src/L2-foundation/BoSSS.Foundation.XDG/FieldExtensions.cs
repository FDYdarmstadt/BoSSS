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


namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Extension methods for <see cref="BoSSS.Foundation.DGField"/>
    /// </summary>
    public static class FieldExtensions {

        /// <summary>
        /// accumulates (to field <paramref name="f"/>), in every cell the
        /// distance form the cut cells times <paramref name="alpha"/>;
        /// Note: this is NOT the geometric distance, but the distance index
        /// that identifies the 'Near+1' , 'Near-1', 'Near+2',  ..., - layers.
        /// </summary>
        public static void AccLevelSetDist(this DGField f, double alpha, LevelSetTracker LevSetTrk, int LevSetIdx) {
            int J = f.Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

            for (int j = 0; j < J; j++) {
                int dist = LevelSetTracker.DecodeLevelSetDist(LevSetTrk._Regions.m_LevSetRegions[j], 0);
                f.SetMeanValue(j, f.GetMeanValue(j) + dist * alpha);
            }
        }

        /// <summary>
        /// accumulates (to field <paramref name="f"/>), in every cell the
        /// number of species times <paramref name="alpha"/>
        /// </summary>
        public static void AccNoOfSpecies(this DGField f, double alpha, LevelSetTracker LevSetTrk, int LevSetIdx) {
            int J = f.Basis.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

            for (int j = 0; j < J; j++) {
                ReducedRegionCode rrc;
                int NoOfSpec = LevSetTrk.GetNoOfSpecies(j, out rrc);
                f.SetMeanValue(j, f.GetMeanValue(j) + NoOfSpec * alpha);
            }
        }
    }
}
