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

namespace BoSSS.Solution.Statistic.QuadRules {

    /// <summary>
    /// Interface for shape classes that represent the surface.
    /// Requirement for the generation of nodes and weights
    /// for integration along the surface
    /// </summary>
    public interface ISurfaceEvaluation {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="testnodes"></param>
        /// <param name="quadwghts"></param>
        /// <param name="nw"></param>
        void CreateNodesAndWeights(out double[,] testnodes, out double[] quadwghts, int nw);

    }

}
