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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.CompressibleFlowCommon.Boundary {

    /// <summary>
    /// Boundary condition mapping for compressible XDG methods
    /// </summary>
    public class XDGCompressibleBoundaryCondMap : CompressibleBoundaryCondMap {

        private string[] speciesNames;

        public Dictionary<string, BoundaryCondition>[] XDGConditionMap {
            get;
            private set;
        }

        public XDGCompressibleBoundaryCondMap(IGridData gridData, AppControl control, Material material, string[] speciesNames) : base(gridData, control, material) {
            this.speciesNames = speciesNames;
            InitXDGConditionMap();
        }

        private void InitXDGConditionMap() {
            // Init
            XDGConditionMap = new Dictionary<string, BoundaryCondition>[GridCommons.FIRST_PERIODIC_BC_TAG];

            foreach (byte edgeTag in gridData.EdgeTagNames.Keys) {
                // Only process non-periodic boundary edges
                if (edgeTag == 0 || edgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                    continue;
                }

                if (XDGConditionMap[edgeTag] != null) {
                    // Already dealt with
                    continue;
                }

                CompressibleBcType bcType = base.EdgeTag2Type[edgeTag];

                //if (!(bcType == CompressibleBcType.adiabaticSlipWall)
                //    && !(bcType == CompressibleBcType.adiabaticWall)
                //    && !(bcType == CompressibleBcType.isothermalWall)
                //    && !(bcType == CompressibleBcType.ringleb)
                //    && !(bcType == CompressibleBcType.symmetryPlane)
                //    ) 
                {
                    foreach (string speciesName in this.speciesNames) {
                        // Init if not already initialized
                        if (XDGConditionMap[edgeTag] == null) {
                            XDGConditionMap[edgeTag] = new Dictionary<string, BoundaryCondition>();
                        }
                        XDGConditionMap[edgeTag].Add(speciesName, BoundaryConditionFactory(bcType, Material, edgeTag, speciesName));
                    }
                }
            }
        }

        public StateVector GetBoundaryState(byte EdgeTag, double time, double[] x, double[] normal, StateVector stateVector, string speciesName) {
            return XDGConditionMap[EdgeTag][speciesName].GetBoundaryState(time, x, normal, stateVector);
        }
    }


}
