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
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;
using System.Diagnostics;
using ilPSP;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Foundation;
using IntersectingQuadrature.TensorAnalysis;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.CompressibleFlowCommon.Boundary {

    /// <summary>
    /// Implementation of the boundary condition for a supersonic inlet (i.e.
    /// Ma > 1.0). In this case, all values at the boundary have to be
    /// prescribed. For this matter, we choose the density, the momentum and
    /// the pressure as input.
    /// </summary>
    public class SpaceTimeBoundaryCondition : BoundaryCondition {

        /// <summary>
        /// field storing the vlaues from previous time step
        /// </summary>
        public DGField[] previous_U;
        
        public SpaceTimeBoundaryCondition(MaterialProperty.Material config)
            : base(config) {
        }

        /// <summary>
        /// Calculates the boundary values for a supersonic inlet (Mach number
        /// greater than 1). Theoretically, we have to impose five conditions
        /// (in the inviscid as well as in the viscid case!) which is why
        /// calculate the complete state from the given density, momentum and
        /// pressure
        /// </summary>
        /// <returns>
        /// \f$ (\rho^*, u_0^*[, u_1^*[, u_2^*]], p*)^T\f$ 
        /// </returns>
        public override StateVector GetBoundaryState(double time, Vector x, Vector normal, StateVector stateIn) {
            Console.WriteLine(" SpaceTimeBoundaryCondtion: GetBoundaryState:  This should never be executed");
            return stateIn;
        }



        /// <summary>
        /// Vectorized implementation of <see cref="GetBoundaryState(double, Vector, Vector, StateVector)"/>
        /// </summary>
        public override void GetBoundaryState(MultidimensionalArray[] StateOut, double time, MultidimensionalArray X, MultidimensionalArray Normals, MultidimensionalArray[] StateIn, int Offset, int NoOfEdges, bool normalFlipped, Material material) {
            Console.WriteLine(" SpaceTimeBoundaryCondtion: GetBoundaryState:  This should never be executed");
        }
    }
}
