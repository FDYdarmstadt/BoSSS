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

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using ilPSP;
using System;
using System.Linq;

namespace CNS.IBM {

    /// <summary>
    /// Variant of <see cref="Convection.ConvectiveCFLConstraint"/> for flows with
    /// immersed boundaries
    /// </summary>
    public class IBMConvectiveCFLConstraint : CFLConstraint {

        private CNSControl config;

        private ImmersedSpeciesMap speciesMap;

        /// <summary>
        /// Constructs a constraint that respects the given
        /// <paramref name="speciesMap"/>
        /// </summary>
        /// <param name="config"></param>
        /// <param name="gridData"></param>
        /// <param name="workingSet"></param>
        /// <param name="speciesMap"></param>
        public IBMConvectiveCFLConstraint(
            CNSControl config, GridData gridData, CNSFieldSet workingSet, ISpeciesMap speciesMap)
            : base(gridData, workingSet) {

            this.config = config;
            this.speciesMap = speciesMap as ImmersedSpeciesMap;

            if (speciesMap == null) {
                throw new ArgumentException(
                    "This type requires an instance of 'ImmersedSpeciesMap'",
                    "speciesMap");
            }
        }

        /// <summary>
        /// Mostly identical to
        /// <see cref="Convection.ConvectiveCFLConstraint.GetCFLStepSize"/>,
        /// but uses the level set to omit nodes in the void domain.
        /// </summary>
        /// <param name="i0"></param>
        /// <param name="Length"></param>
        /// <returns>
        /// <see cref="Convection.ConvectiveCFLConstraint.GetCFLStepSize"/>
        /// </returns>
        protected override double GetCFLStepSize(int i0, int Length) {
            int iKref = gridData.Cells.GetRefElementIndex(i0);
            int noOfNodesPerCell = base.EvaluationPoints[iKref].NoOfNodes;
            int D = gridData.Grid.SpatialDimension;
            Material material = speciesMap.GetMaterial(double.NaN);

            MultidimensionalArray levelSetValues =
                speciesMap.Tracker.DataHistories[0].Current.GetLevSetValues(base.EvaluationPoints[iKref], i0, Length);

            SpeciesId species = speciesMap.Tracker.GetSpeciesId(speciesMap.Control.FluidSpeciesName);
            //var hMinArray = speciesMap.QuadSchemeHelper.CellAgglomeration.CellLengthScales[species];
            var volFrac = speciesMap.CellAgglomeration.CellVolumeFrac[species];
            var hMin = gridData.Cells.h_min;

            double cfl = double.MaxValue;
            for (int i = 0; i < Length; i++) {
                int cell = i0 + i;

                // Option 1: Use volume fraction times traditional metric, such
                // that time-steps are identical to non-IBM cases when no
                // interface is present
                double hmin = hMin[cell] * volFrac[cell];

                // Option 2: Use length scale "volume over surface", which seems
                // to be more robust for awkward cuts. However, yields
                // significantly smaller time-steps in uncut cells
                //double hmin = hMinArray[cell];


                for (int node = 0; node < noOfNodesPerCell; node++) {
                    if (levelSetValues[i, node].Sign() != (double)speciesMap.Control.FluidSpeciesSign) {
                        continue;
                    }

                    Vector3D momentum = new Vector3D();
                    for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                        momentum[d] = momentumValues[d][i, node];
                    }



                    StateVector state = new StateVector(
                        material, densityValues[i, node], momentum, energyValues[i, node]);
                    double cflhere = hmin / (state.Velocity.Abs() + state.SpeedOfSound);

#if DEBUG
                    if (double.IsNaN(cflhere)) {
                        throw new Exception("Could not determine CFL number");
                    }
#endif

                    cfl = Math.Min(cfl, cflhere);
                }
            }

            int degree = workingSet.ConservativeVariables.Max(f => f.Basis.Degree);
            int twoNPlusOne = 2 * degree + 1;
            return cfl / twoNPlusOne;
        }
    }
}