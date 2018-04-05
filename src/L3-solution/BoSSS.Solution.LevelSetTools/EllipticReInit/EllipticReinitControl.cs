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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools.Reinit;
using ilPSP.LinSolvers;

namespace BoSSS.Solution.LevelSetTools.EllipticReInit {
    /// <summary>
    /// Implement this in a Control Class, when using the Elliptic ReInit
    /// </summary>
    [Serializable]
    public class EllipticReInitAlgoControl  {

        /// <summary>
        /// L2-Criteria for Change Rate
        /// </summary>
        public double ConvergenceCriterion = 1e-4;
        /// <summary>
        /// Maximal Number of Iterations
        /// </summary>
        public int MaxIt = 1000;
        
        /// <summary>
        /// Potential Function for ReInit
        /// DoubleWell by Basting Seems to work best so far
        /// </summary>
        public ReInitPotential Potential = ReInitPotential.BastingDoubleWell;

        ///// <summary>
        ///// The Variant for Quadrature on the Cut-Cells
        ///// </summary>
        //public XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Classic;

        /// <summary>
        /// Use a fast marching solver to generate an initial solution
        /// </summary>
        public bool FastMarchingPrecond = false;

        /// <summary>
        /// Underrelaxation-Factor for the fix-point iteration
        /// </summary>
        public double underrelaxation=1;

        /// <summary>
        /// Upwind Based Flux Formulation
        /// </summary>
        public bool Upwinding = false;


        /// <summary>
        /// Generates the linear Solver;
        /// </summary>
        public Func<ISparseSolver> solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();


        public double PenaltyMultiplierFlux = 1;

        public double PenaltyMultiplierInterface = 1;

        /// <summary>
        /// Verbose Output - prints Changerate of ReInit-Iterations
        /// </summary>
        public bool PrintIterations = false;

    }


}
