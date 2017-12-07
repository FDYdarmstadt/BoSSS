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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using CNS.EquationSystem;
using CNS.IBM;
using ilPSP;
using System;
using System.Diagnostics;
using System.Linq;

namespace CNS.ShockCapturing {

    /// <summary>
    /// Implementation of the CFL constraint induced by the diffusive terms of
    /// the compressible Navier-Stokes equations. The exact implementation is
    /// based on the formulas presented in Gassner, Loercher, Munz 2008:
    /// A Discontinuous Galerkin Scheme based on a Space-Time Expansion II.
    /// Viscous Flow Equations in Multi Dimensions
    /// </summary>
    public class ArtificialViscosityCFLConstraint : CFLConstraint {

        public ArtificialViscosityCFLConstraint(
            CNSControl config, GridData gridData, CNSFieldSet workingSet, ISpeciesMap speciesMap)
            : base(gridData, workingSet) {

            this.config = config;
            this.speciesMap = speciesMap;

            if (gridData.Grid.RefElements.Length > 1) {
                throw new NotImplementedException();
            }
        }

        private CNSControl config;

        private ISpeciesMap speciesMap;

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
            double scaling = Math.Max(4.0 / 3.0, config.EquationOfState.HeatCapacityRatio / config.PrandtlNumber);
            DGField artificialViscosity = workingSet.ParameterFields.Where(c => c.Identification.Equals(Variables.ArtificialViscosity)).Single();
            var hmin = gridData.Cells.h_min;

            double cfl = double.MaxValue;
            switch (speciesMap) {
                case ImmersedSpeciesMap ibmMap: {
                        MultidimensionalArray levelSetValues = ibmMap.Tracker.GetLevSetValues(0, base.EvaluationPoints[iKref], i0, Length);
                        SpeciesId species = ibmMap.Tracker.GetSpeciesId(ibmMap.Control.FluidSpeciesName);
                        var hMinArray = ibmMap.QuadSchemeHelper.CellAgglomeration.CellLengthScales[species];
                        
                        for (int i = 0; i < Length; i++) {
                            int cell = i0 + i;

                            for (int node = 0; node < noOfNodesPerCell; node++) {
                                if (levelSetValues[i, node].Sign() != (double)ibmMap.Control.FluidSpeciesSign) {
                                    continue;
                                }
                                
                                double hminlocal = hMinArray[cell];

                                double nu = artificialViscosity.GetMeanValue(cell) / config.ReynoldsNumber;
                                Debug.Assert(!double.IsNaN(nu), "IBMArtificialViscosityCFLConstraint: nu is NaN!");

                                double cflhere;
                                if (nu == 0) {
                                    cflhere = double.MaxValue;
                                } else {
                                    cflhere = hminlocal * hminlocal / scaling / nu;
                                }

#if DEBUG
                                if (double.IsNaN(cflhere)) {
                                    throw new Exception("Could not determine CFL number");
                                }
#endif

                                cfl = Math.Min(cfl, cflhere);
                            }
                        }
                    }
                    break;

                default: {
                        for (int i = 0; i < Length; i++) {
                            int cell = i0 + i;

                            for (int node = 0; node < noOfNodesPerCell; node++) {
                                double hminlocal = gridData.Cells.h_min[cell];

                                double nu = artificialViscosity.GetMeanValue(cell) / config.ReynoldsNumber;
                                Debug.Assert(!double.IsNaN(nu), "ArtificialViscosityCFLConstraint: nu is NaN!");

                                double cflhere;
                                if (nu == 0) {
                                    cflhere = double.MaxValue;
                                } else {
                                    cflhere = hminlocal * hminlocal / scaling / nu;
                                }

                                cfl = Math.Min(cfl, cflhere);
                            }
                        }
                    }
                    break;
            }

            int degree = workingSet.ConservativeVariables.Max(f => f.Basis.Degree);
            int twoNPlusOne = 2 * degree + 1;
            return cfl * GetBetaMax(degree) / twoNPlusOne / twoNPlusOne / Math.Sqrt(CNSEnvironment.NumberOfDimensions);
        }
    }
}
