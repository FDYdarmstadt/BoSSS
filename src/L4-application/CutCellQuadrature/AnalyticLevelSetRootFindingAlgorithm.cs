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
using BoSSS.Foundation.XDG.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.Quadrature.HMF;

namespace CutCellQuadrature {

    class AnalyticLevelSetRootFindingAlgorithm : LineSegment.IRootFindingAlgorithm {

        #region IRootFindingAlgorithm Members

        public double Tolerance {
            get {
                return 1e-16;
            }
        }

        public double[] GetRoots(LineSegment segment, ILevelSet levelSet, int cell, int iKref) {
            if (levelSet is IAnalyticLevelSet) {
                return ((IAnalyticLevelSet)levelSet).GetRoots(segment, cell).ToArray();
            } else {
                throw new ArgumentException("Level set must implement IAnalyticLevelSet", "levelSet");
            }
        }

        #endregion
    }
}
