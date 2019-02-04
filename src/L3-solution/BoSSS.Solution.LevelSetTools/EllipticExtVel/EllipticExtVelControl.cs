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
using BoSSS.Solution;
using BoSSS.Foundation.XDG;
using static BoSSS.Solution.AdvancedSolvers.DirectSolver;
using ilPSP.LinSolvers;

namespace BoSSS.Solution.LevelSetTools.EllipticExtension {

    
    public enum FluxVariant {
        GradientBased,
        ValueBased,
        SWIP
    }

    /// <summary>
    /// Implement this in a Control Class, when using the Elliptic ReInit
    /// </summary>
    [System.Serializable]
    public class EllipticExtVelAlgoControl {

        
        /// <summary>
        /// Switch for the choice of the flux directions
        /// </summary>
        public FluxVariant FluxVariant = FluxVariant.GradientBased;

        public bool subGridRestriction = false;

        /// <summary>
        /// Mutliplier for the penalty factor  between the cells.
        /// </summary>
        public double PenaltyMultiplierFlux = 1;

        /// <summary>
        /// Mutliplier for the penalty factor at the interface
        /// </summary>
        public double PenaltyMultiplierInterface = 1;

        /// <summary>
        /// LinearSolver for the Extension Problem
        /// </summary>
        public Func<ISparseSolver> solverFactory = () => new ilPSP.LinSolvers.MUMPS.MUMPSSolver();


        /// <summary>
        /// The prefactor for isotropic viscosity - this should be chosen as small as possible
        /// </summary>
        public double IsotropicViscosity = 0;

        //using (var slv = new ilPSP.LinSolvers.HYPRE.GMRES()) {
        //using (var slv = new ilPSP.LinSolvers.MUMPS.MUMPSSolver()) {
        //slv.Tolerance = 1E-10;
        //slv.MaxIterations = 10000;
        //slv.PrintLevel = 2;
        //slv.KrylovSpaceMaxDimension = 20;
    }
}
