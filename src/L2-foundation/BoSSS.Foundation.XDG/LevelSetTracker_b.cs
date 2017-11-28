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
using System.Collections;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP.Utils;
using System.Diagnostics;
using System.Linq;
using ilPSP;
using BoSSS.Foundation.Comm;
using System.Collections.Generic;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.XDG {

    partial class LevelSetTracker {

        

        public Dictionary<XQuadFactoryHelper.MomentFittingVariants, XQuadFactoryHelper> m_QuadFactoryHelpers
            = new Dictionary<XQuadFactoryHelper.MomentFittingVariants,XQuadFactoryHelper>();

        /*
        /// <summary>
        /// Central 'factory' for creating Level Set - related quadrature.
        /// </summary>
        /// <remarks>
        /// The centralized approach should avoid multiple creation of the same quadrature rule.
        /// </remarks>
        public XQuadFactoryHelper GetXQuadFactoryHelper(XQuadFactoryHelper.MomentFittingVariants variant) {
            if (!m_QuadFactoryHelpers.ContainsKey(variant)) {
                m_QuadFactoryHelpers[variant] = new XQuadFactoryHelper(this, variant);
            }
            return m_QuadFactoryHelpers[variant];
        }
        
      
        */


        public XDGSpaceMetrics GetXDGSpaceMetrics(XQuadFactoryHelper.MomentFittingVariants variant, int quadorder, int stackindex) {

        }
    }
}
