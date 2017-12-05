using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using CNS.EquationSystem;
using CNS.MaterialProperty;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CNS.IBM {
    /// <summary>
    /// Implementation of the CFL constraint induced by the diffusive terms of
    /// the compressible Navier-Stokes equations. The exact implementation is
    /// based on the formulas presented in Gassner, Loercher, Munz 2008:
    /// A Discontinuous Galerkin Scheme based on a Space-Time Expansion II.
    /// Viscous Flow Equations in Multi Dimensions.
    /// <seealso cref="IBMDiffusiveCFLConstraint"/>
    /// </summary>
    public class IBMArtificialViscosityCFLConstraint : CFLConstraint {

        /// <summary>
        /// Factors by Gassner, see <see cref="GetBetaMax"/>
        /// </summary>
        private static double[] beta_max = new double[] {
            // 2.0 is just a guess
            2.0, 1.46, 0.8, 0.54, 0.355, 0.28, 0.21, 0.16, 0.14, 0.12, 0.1 };

        private CNSControl config;

        private ImmersedSpeciesMap speciesMap;

        /// <summary>
        /// Constructs a new CFL constraint based on artificial viscosity
        /// for IBM simulations
        /// </summary>
        /// <param name="config"></param>
        /// <param name="gridData"></param>
        /// <param name="workingSet"></param>
        /// <param name="speciesMap"></param>
        public IBMArtificialViscosityCFLConstraint(
            CNSControl config, GridData gridData, CNSFieldSet workingSet, ISpeciesMap speciesMap)
            : base(gridData, workingSet) {

            this.config = config;
            this.speciesMap = speciesMap as ImmersedSpeciesMap;

            if (gridData.Grid.RefElements.Length > 1) {
                throw new NotImplementedException();
            }
        }

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

            // IBM part
            MultidimensionalArray levelSetValues = speciesMap.Tracker.GetLevSetValues(0, base.EvaluationPoints[iKref], i0, Length);
            SpeciesId species = speciesMap.Tracker.GetSpeciesId(speciesMap.Control.FluidSpeciesName);
            var hMinArray = speciesMap.QuadSchemeHelper.CellAgglomeration.CellLengthScales[species];

            double cfl = double.MaxValue;
            for (int i = 0; i < Length; i++) {
                int cell = i0 + i;

                for (int node = 0; node < noOfNodesPerCell; node++) {
                    if (levelSetValues[i, node].Sign() != (double)speciesMap.Control.FluidSpeciesSign) {
                        continue;
                    }

                    // IBM part
                    double hmin = hMinArray[cell];

                    DGField artificialViscosity = workingSet.ParameterFields.Where(c => c.Identification.Equals("artificialViscosity")).Single();
                    double nu = artificialViscosity.GetMeanValue(cell) / config.ReynoldsNumber;
                    Debug.Assert(!double.IsNaN(nu), "IBMArtificialViscosityCFLConstraint: nu is NaN!");

                    double coeff = Math.Max(4.0 / 3.0, config.EquationOfState.HeatCapacityRatio / config.PrandtlNumber);

                    double cflhere;
                    if (nu == 0) {
                        cflhere = double.MaxValue;
                    } else {
                        cflhere = hmin * hmin / coeff / nu;
                    }

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
            return cfl * GetBetaMax(degree) / twoNPlusOne / twoNPlusOne / Math.Sqrt(CNSEnvironment.NumberOfDimensions);
        }
    }
}
