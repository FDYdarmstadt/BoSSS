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

using BoSSS.Foundation.XDG;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Collection of various extension methods.
    /// </summary>
    static public class ExtensionMethods {

        /// <summary>
        /// Updates the XDG component of an aggregation basis according to the new agglomeration.
        /// </summary>
        public static void UpdateXdgAggregationBasis(this AggregationGridBasis[][] MultigridBasis, MultiphaseCellAgglomerator CurrentAgglomeration) {

            bool useX = false;
            var _XdgAggregationBasis = new XdgAggregationBasis[MultigridBasis.Length];
            for (int iLevel = 0; iLevel < MultigridBasis.Length; iLevel++) {
                XdgAggregationBasis[] xab  = MultigridBasis[iLevel].Where(b => b is XdgAggregationBasis).Select(b => ((XdgAggregationBasis)b)).ToArray();
                if (xab != null && xab.Length > 0) {
                    for (int ib = 1; ib < xab.Length; ib++) {
                        if (!(object.ReferenceEquals(xab[ib].DGBasis.GridDat, CurrentAgglomeration.Tracker.GridDat)))
                            throw new ApplicationException();
                        if (!object.ReferenceEquals(xab[0], xab[ib]))
                            throw new ArgumentException("One should only use one XDG aggregation basis per multigrid level.");
                    }
                    _XdgAggregationBasis[iLevel] = xab[0];
                    useX = true;
                }
            }

            if (useX) {
                foreach (var xmgb in _XdgAggregationBasis) {
                    xmgb.Update(CurrentAgglomeration);
                }
            }
        }
    }
}
