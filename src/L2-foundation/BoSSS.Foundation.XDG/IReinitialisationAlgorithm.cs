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
namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// common template of all Level-Set reinitialization algorithms
    /// </summary>
    public interface IReinitialisationAlgorithm {


        /// <summary>
        /// performs a reinitialization of the Level Set function <paramref name="LS"/>
        /// </summary>
        /// <param name="LS">
        /// on entry, an arbitrary Level Set function; on exit, a signed-distance Level Set function 
        /// with equal zero-set.
        /// </param>
        /// <param name="Restriction">
        /// optional restriction to the computational domain
        /// </param>
        void ReInitialize(LevelSet LS, SubGrid Restriction=null);

    }
}
