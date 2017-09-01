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
using System.Text;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.Quadrature;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    /// <summary>
    /// Common interface for analytically prescribed level sets
    /// </summary>
    public interface IAnalyticLevelSet : ILevelSet {

        /// <summary>
        /// Determines the roots of the level set on
        /// <paramref name="lineSegment"/> in the given
        /// <paramref name="element"/>.
        /// </summary>
        /// <param name="lineSegment">
        /// The considered line segment
        /// </param>
        /// <param name="element">
        /// The index of the considered element
        /// </param>
        /// <returns>
        /// The parameter values t of the parametrization of
        /// <paramref name="lineSegment"/> for which this level set is zero in
        /// element <paramref name="element"/>
        /// </returns>
        IEnumerable<double> GetRoots(LineSegment lineSegment, int element);
    }
}
