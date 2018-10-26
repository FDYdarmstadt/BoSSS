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

using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;

namespace CNS.Boundary {

    /// <summary>
    /// Base class for implementations of boundary value calculations. That is,
    /// subclasses of this class determine the characteristic values at the
    /// outside of the flow domain (i.e. the "+"-values for edge fluxes) which
    /// may then be used to weakly impose boundary conditions
    /// </summary>
    abstract public class BoundaryCondition {

        /// <summary>
        /// Configuration options
        /// </summary>
        protected CNSControl config;

        /// <summary>
        /// Constructs a new boundary condition
        /// </summary>
        /// <param name="config">Configuration options</param>
        protected BoundaryCondition(CNSControl config) {
            this.config = config;
        }

        /// <summary>
        /// Implement this method in order to calculate the values of the flow
        /// variables (<see cref="StateVector"/>) outside the computational
        /// domain (i.e. on "+"-side of a boundary edge).
        /// </summary>
        /// <param name="time">The solution time</param>
        /// <param name="x">The position in global coordinates</param>
        /// <param name="normal">The outward normal vector</param>
        /// <param name="stateIn">The values on the inside of the cell</param>
        /// <returns>The values on the outside of the cell</returns>
        abstract public StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn);

        /// <summary>
        /// Utility function create a <see cref="Vector"/> out of double[]
        /// representation of the outward unit normal.
        /// </summary>
        /// <param name="normal">The original normal to be transformed</param>
        /// <returns>
        /// A 3D vector representation of the normal. If
        /// <paramref name="normal"/>.Length is smaller than three the
        /// remaining elements of the vector are set to zero.
        /// </returns>
        protected static Vector GetNormalVector(double[] normal) {
            Vector normalVector = new Vector();
            for (int i = 0; i < normal.Length; i++) {
                normalVector[i] = normal[i];
            }
            return normalVector;
        }
    }
}
