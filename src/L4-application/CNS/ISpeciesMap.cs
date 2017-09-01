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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using CNS.MaterialProperty;

namespace CNS {

    /// <summary>
    /// A map that determines the active species in some point.
    /// </summary>
    public interface ISpeciesMap {

        /// <summary>
        /// Determines the material/fluid (see
        /// <see cref="IEquationOfState"/>) of the active species in some point
        /// of the domain depending on the value of level set.
        /// </summary>
        /// <param name="levelSetValue">
        /// The value of the level set function that determines th active
        /// species.
        /// </param>
        /// <returns>
        /// The valid equation of state for the given
        /// <paramref name="levelSetValue"/>.
        /// </returns>
        Material GetMaterial(double levelSetValue);

        /// <summary>
        /// Returns the sub-domain associated with the represented species
        /// </summary>
        SubGrid SubGrid {
            get;
        }

        /// <summary>
        /// Information about the grid
        /// </summary>
        GridData GridData {
            get;
        }
    }
}
