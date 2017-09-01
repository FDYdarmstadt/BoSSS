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
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using CNS.EquationSystem;
using CNS.Exception;
using CNS.MaterialProperty;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace CNS.IBM {

    /// <summary>
    /// Variant of <see cref="Convection.ConvectiveCFLConstraint"/> for flows with
    /// immersed boundaries
    /// </summary>
    public class IBMConvectiveCFLConstraint : CFLConstraint {

        private CNSControl config;

        private ImmersedSpeciesMap speciesMap;

        /// <summary>
        /// Factors by Gassner
        /// </summary>
        private static double[] alpha_max = new double[] {
            // 2.0 is just a guess
            1.8, 1.3, 1.1, 0.9, 0.7, 0.7};

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
                speciesMap.Tracker.GetLevSetValues(0, base.EvaluationPoints[iKref], i0, Length);

            SpeciesId species = speciesMap.Tracker.GetSpeciesId(speciesMap.Control.FluidSpeciesName);
            var hMinArray = speciesMap.QuadSchemeHelper.CellAgglomeration.CellLengthScales[species];
            //var volFrac = speciesMap.QuadSchemeHelper.CellAgglomeration.CellVolumeFrac[species];
            //var hMinGass = speciesMap.h_min;
            //var hMin = gridData.Cells.h_min;

            double cfl = double.MaxValue;
            for (int i = 0; i < Length; i++) {
                int cell = i0 + i;

                //double hmin = hMin[cell] * volFrac[cell];
                double hmin = hMinArray[cell];
                //double hmin = hMinGass[cell];


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
                        throw new NumericalAlgorithmException("Could not determine CFL number");
                    }
#endif

                    cfl = Math.Min(cfl, cflhere);
                }
            }

            int degree = workingSet.ConservativeVariables.Max(f => f.Basis.Degree);
            int twoNPlusOne = 2 * degree + 1;
            //return cfl * alpha_max[degree] / twoNPlusOne;
            return cfl / twoNPlusOne;
        }
    }
}