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
using BoSSS.Platform.LinAlg;
using CNS.MaterialProperty;
using CNS.EquationSystem;
using BoSSS.Foundation.Grid.Classic;

namespace CNS.Convection {

    /// <summary>
    /// Implementation of the CFL constraint induced by the convective terms of
    /// the compressible Navier-Stokes equations
    /// </summary>
    public class ConvectiveCFLConstraint : CFLConstraint {

        /// <summary>
        /// <see cref="ConvectiveCFLConstraint.ConvectiveCFLConstraint"/>
        /// </summary>
        private CNSControl config;

        /// <summary>
        /// <see cref="ConvectiveCFLConstraint.ConvectiveCFLConstraint"/>
        /// </summary>
        private ISpeciesMap speciesMap;

        /// <summary>
        /// Constructs a new constraint
        /// </summary>
        /// <param name="config"></param>
        /// <param name="gridData"></param>
        /// <param name="workingSet"></param>
        /// <param name="speciesMap"></param>
        public ConvectiveCFLConstraint(
            CNSControl config, GridData gridData, CNSFieldSet workingSet, ISpeciesMap speciesMap)
            : base(gridData, workingSet) {

            this.config = config;
            this.speciesMap = speciesMap;

            if (gridData.Grid.RefElements.Length > 1) {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Computes the maximum admissible step-size according to
        /// \f$ 
        /// \text{cfl} = \frac{1}{2N+1} \min \frac{h}{|\vec{u}| + a} 
        /// \f$ 
        /// where \f$ h\f$  is the minimum cell diameter,
        /// \f$ \vec{u}\f$  is the nodal velocity vector
        /// and \f$ a\f$  is the nodal speed of sound.
        /// Additionally, \f$ N\f$  is the degree of the
        /// DG with the highest degree of approximation.
        /// </summary>
        /// <param name="i0"></param>
        /// <param name="Length"></param>
        /// <returns></returns>
        protected override double GetCFLStepSize(int i0, int Length) {
            int iKref = gridData.Cells.GetRefElementIndex(i0);
            int noOfNodesPerCell = EvaluationPoints[iKref].NoOfNodes;
            Material material = speciesMap.GetMaterial(double.NaN);

            double cfl = double.MaxValue;
            for (int i = 0; i < Length; i++) {
                int cell = i0 + i;

                for (int node = 0; node < noOfNodesPerCell; node++) {
                    Vector3D momentum = new Vector3D();
                    for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                        momentum[d] = momentumValues[d][i, node];
                    }

                    StateVector state = new StateVector(
                        material, densityValues[i, node], momentum, energyValues[i, node]);
                    double cflhere = gridData.Cells.h_min[cell] /
                        (state.Velocity.Abs() + material.EquationOfState.GetSpeedOfSound(state));

                    cfl = Math.Min(cfl, cflhere);
                }
            }

            int twoNPlusOne = 2 * workingSet.ConservativeVariables.Max(f => f.Basis.Degree) + 1;
            return cfl / twoNPlusOne;
        }
    }
}