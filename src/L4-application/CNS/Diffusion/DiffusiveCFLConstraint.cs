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
using CNS.IBM;
using CNS.MaterialProperty;
using ilPSP;
using System;
using System.Diagnostics;
using System.Linq;

namespace CNS.Diffusion {

    /// <summary>
    /// Implementation of the CFL constraint induced by the diffusive terms of
    /// the compressible Navier-Stokes equations. The exact implementation is
    /// based on the formulas presented in Gassner, Loercher, Munz 2008:
    /// A Discontinuous Galerkin Scheme based on a Space-Time Expansion II.
    /// Viscous Flow Equations in Multi Dimensions
    /// </summary>
    public class DiffusiveCFLConstraint : CFLConstraint {

        /// <summary>
        /// Factors by Gassner, see <see cref="GetBetaMax"/>
        /// </summary>
        private static double[] beta_max = new double[] {
            // 2.0 is just a guess
            2.0, 1.46, 0.8, 0.54, 0.355, 0.28, 0.21, 0.16, 0.14, 0.12, 0.1, 0.08 };

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

        private ISpeciesMap speciesMap;

        /// <summary>
        /// Constructs a new constraint
        /// </summary>
        /// <param name="config"></param>
        /// <param name="gridData"></param>
        /// <param name="workingSet"></param>
        /// <param name="speciesMap"></param>
        public DiffusiveCFLConstraint(
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
        /// GassnerEtAl2008, equation 67.
        /// </summary>
        /// <param name="i0"></param>
        /// <param name="Length"></param>
        /// <returns></returns>
        protected override double GetCFLStepSize(int i0, int Length) {
            int iKref = gridData.Cells.GetRefElementIndex(i0);
            int noOfNodesPerCell = base.EvaluationPoints[iKref].NoOfNodes;
            Material material = speciesMap.GetMaterial(double.NaN);
            var hmin = gridData.Cells.h_min;
            double scaling = Math.Max(4.0 / 3.0, config.EquationOfState.HeatCapacityRatio / config.PrandtlNumber);

            // Following code is performance-critical, so expect spagetti code ahead
            double cfl = double.MaxValue;
            switch (speciesMap) {
                // 1) IBM optimized (ignores nodes in the void part of a cell, only applies to constant viscosity gases)
                case ImmersedSpeciesMap ibmMap when material.ViscosityLaw is ConstantViscosity: {
                        MultidimensionalArray levelSetValues =
                            ibmMap.Tracker.GetLevSetValues(0, base.EvaluationPoints[iKref], i0, Length);

                        SpeciesId species = ibmMap.Tracker.GetSpeciesId(ibmMap.Control.FluidSpeciesName);
                        //var hMin = speciesMap.QuadSchemeHelper.CellAgglomeration.CellLengthScales[species];
                        var volFrac = ibmMap.QuadSchemeHelper.CellAgglomeration.CellVolumeFrac[species];

                        double nu = 1.0 / config.ReynoldsNumber;
                        for (int i = 0; i < Length; i++) {
                            int cell = i0 + i;

                            double hminlocal = hmin[cell] * volFrac[cell];
                            //double hminlocal = hMin[cell];
                            //double hmin = hMinGass[cell];

                            for (int node = 0; node < noOfNodesPerCell; node++) {
                                if (levelSetValues[i, node].Sign() != (double)ibmMap.Control.FluidSpeciesSign) {
                                    continue;
                                }

                                double cflhere = hminlocal * hminlocal / scaling / nu;

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

                // 2) IBM generic (ignores nodes in the void part of a cell)
                case ImmersedSpeciesMap ibmMap: {
                        MultidimensionalArray levelSetValues =
                            ibmMap.Tracker.GetLevSetValues(0, base.EvaluationPoints[iKref], i0, Length);

                        SpeciesId species = ibmMap.Tracker.GetSpeciesId(ibmMap.Control.FluidSpeciesName);
                        //var hMin = speciesMap.QuadSchemeHelper.CellAgglomeration.CellLengthScales[species];
                        var volFrac = ibmMap.QuadSchemeHelper.CellAgglomeration.CellVolumeFrac[species];

                        for (int i = 0; i < Length; i++) {
                            int cell = i0 + i;

                            double hminlocal = hmin[cell] * volFrac[cell];
                            //double hminlocal = hMin[cell];
                            //double hmin = hMinGass[cell];

                            for (int node = 0; node < noOfNodesPerCell; node++) {
                                if (levelSetValues[i, node].Sign() != (double)ibmMap.Control.FluidSpeciesSign) {
                                    continue;
                                }

                                Vector3D momentum = new Vector3D();
                                for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                                    momentum[d] = momentumValues[d][i, node];
                                }

                                StateVector state = new StateVector(
                                    material, densityValues[i, node], momentum, energyValues[i, node]);
                                double nu = state.GetViscosity(cell) / config.ReynoldsNumber;
                                Debug.Assert(nu > 0.0);

                                double cflhere = hminlocal * hminlocal / scaling / nu;

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

                // 3) Non-IBM optimized (only applies to constant viscosity gases)
                case var speciesMap when material.ViscosityLaw is ConstantViscosity: {
                        double nu = 1.0 / config.ReynoldsNumber; // Constant viscosity
                        for (int i = 0; i < Length; i++) {
                            int cell = i0 + i;

                            for (int node = 0; node < noOfNodesPerCell; node++) {
                                double hminlocal = gridData.Cells.h_min[cell];
                                double cflhere = hminlocal * hminlocal / scaling / nu;
                                cfl = Math.Min(cfl, cflhere);
                            }
                        }
                    }
                    break;

                // 4) Non-IBM generic (works for all viscosity laws)
                default: {
                        for (int i = 0; i < Length; i++) {
                            int cell = i0 + i;

                            for (int node = 0; node < noOfNodesPerCell; node++) {
                                Vector3D momentum = new Vector3D();
                                for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                                    momentum[d] = momentumValues[d][i, node];
                                }

                                StateVector state = new StateVector(
                                    material, densityValues[i, node], momentum, energyValues[i, node]);

                                double hminlocal = hmin[cell];

                                double nu = state.GetViscosity(cell) / config.ReynoldsNumber;
                                Debug.Assert(nu > 0.0);

                                double cflhere = hminlocal * hminlocal / scaling / nu;
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
