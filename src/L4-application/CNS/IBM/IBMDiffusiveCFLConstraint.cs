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
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using CNS.Exception;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;

namespace CNS.IBM {

    /// <summary>
    /// Variant of <see cref="Diffusion.DiffusiveCFLConstraint"/> for immersed
    /// boundary runs
    /// </summary>
    public class IBMDiffusiveCFLConstraint : CFLConstraint {

        /// <summary>
        /// Factors by Gassner, see <see cref="GetBetaMax"/>
        /// </summary>
        private static double[] beta_max = new double[] {
            // 2.0 is just a guess
            2.0, 1.46, 0.8, 0.54, 0.355, 0.28, 0.21, 0.16, 0.14, 0.12, 0.1 };

        /// <summary>
        /// Get experimentally obtained stability limit by Gassner for
        /// different polynomial degrees.
        /// </summary>
        /// <param name="polydegree">
        /// The polynomial degree of the approximation
        /// </param>
        /// <returns>
        /// An experimental factor for the stability limit
        /// </returns>
        private static double GetBetaMax(int polydegree) {
            return beta_max[polydegree];
        }

        private CNSControl config;

        private ImmersedSpeciesMap speciesMap;

        /// <summary>
        /// Constructs a new constraint
        /// </summary>
        /// <param name="config"></param>
        /// <param name="gridData"></param>
        /// <param name="workingSet"></param>
        /// <param name="speciesMap"></param>
        public IBMDiffusiveCFLConstraint(
            CNSControl config, GridData gridData, CNSFieldSet workingSet, ISpeciesMap speciesMap)
            : base(gridData, workingSet) {

            this.config = config;
            this.speciesMap = speciesMap as ImmersedSpeciesMap;

            if (gridData.Grid.RefElements.Length > 1) {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Computes the maximum admissible step-size according to
        /// GassnerEtAl2008, equation 67.
        /// </summary>
        /// <param name="i0"></param>
        /// <param name="Length"></param>
        /// <returns></returns>
        protected override double GetCFLStepSize(int i0, int Length) {
            int iKref = gridData.Cells.GetRefElementIndex(i0);
            int noOfNodesPerCell = base.EvaluationPoints[iKref].NoOfNodes;

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

                    Material material = speciesMap.GetMaterial(double.NaN);
                    Vector3D momentum = new Vector3D();
                    for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                        momentum[d] = momentumValues[d][i, node];
                    }

                    StateVector state = new StateVector(
                        material, densityValues[i, node], momentum, energyValues[i, node]);

                    double coeff = Math.Max(4.0 / 3.0, config.EquationOfState.HeatCapacityRatio / config.PrandtlNumber);
                    double cflhere = hmin * hmin / coeff / (state.GetViscosity(cell) / config.ReynoldsNumber);

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
            return cfl * GetBetaMax(degree) / twoNPlusOne / twoNPlusOne / Math.Sqrt(CNSEnvironment.NumberOfDimensions);
            //return cfl / twoNPlusOne / twoNPlusOne;
        }
    }
}
