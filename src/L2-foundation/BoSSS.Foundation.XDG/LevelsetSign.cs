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
    /// used to encode the sign of the level set (in non-performance
    /// critical regions) in one cell or edge;
    /// </summary>
    public enum LevelsetSign {

        /// <summary>
        /// level set is negative in the whole cell/edge, i.e.
        /// </summary>
        Negative = -1,

        /// <summary>
        /// both, positive and negative (usually a cut cell/edge or a cell/edge
        /// in the near-region);
        /// </summary>
        Both = 0,

        /// <summary>
        /// level set is positive
        /// </summary>
        Positive = 1
    }
}
