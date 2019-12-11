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
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using ilPSP;
using System;

namespace BoSSS.Solution.CompressibleFlowCommon.Boundary {

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
        protected MaterialProperty.Material config;

        /// <summary>
        /// Constructs a new boundary condition
        /// </summary>
        /// <param name="config">Configuration options</param>
        protected BoundaryCondition(MaterialProperty.Material config) {
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
        abstract public StateVector GetBoundaryState(double time, Vector x, Vector normal, StateVector stateIn);

        /// <summary>
        /// Vectorized version of <see cref="GetBoundaryState(double, Vector, Vector, StateVector)"/>
        /// </summary>
        public virtual void GetBoundaryState(MultidimensionalArray[] StateOut, double time, MultidimensionalArray X, MultidimensionalArray Normals, MultidimensionalArray[] StateIn, int Offset, int NoOfEdges, bool normalFlipped, MaterialProperty.Material material) {
            if (X.Dimension != 3)
                throw new ArgumentException();
            int D = X.GetLength(2);
            int NoOfNodes = X.GetLength(1);
            double sign = normalFlipped ? -1.0 : 1.0;

            if (StateIn.Length != D + 2)
                throw new ArgumentException();
            if (StateOut.Length != D + 2)
                throw new ArgumentException();
            bool is2D = D >= 2;
            bool is3D = D >= 3;
            if (D < 1 || D > 3)
                throw new NotSupportedException();

            var Density = StateOut[0];
            var Energy = StateOut[D + 1];
            var MomentumX = StateOut[1];
            var MomentumY = is2D ? StateOut[2] : null;
            var MomentumZ = is3D ? StateOut[3] : null;

            Vector xLocal = new Vector(D);
            Vector normalLocal = new Vector(D);
            for (int e = 0; e < NoOfEdges; e++) {
                int edge = e + Offset;

                // Loop over nodes
                for (int n = 0; n < NoOfNodes; n++) {
                    xLocal.x = X[edge, n, 0];
                    normalLocal.x = Normals[edge, n, 0] * sign;
                    if (is2D) {
                        xLocal.y = X[edge, n, 1];
                        normalLocal.y = Normals[edge, n, 1] * sign;
                    }
                    if(is3D) {
                        xLocal.z = X[edge, n, 2];
                        normalLocal.z = Normals[edge, n, 2] * sign;
                    }

                    StateVector stateIn = new StateVector(material, StateIn, edge, n, D);
                    //OptimizedHLLCFlux.State.Start();
                    StateVector stateBoundary = GetBoundaryState(time, xLocal, normalLocal, stateIn);
                    //OptimizedHLLCFlux.State.Stop();

                    Density[edge, n] = stateBoundary.Density;
                    MomentumX[edge, n] = stateBoundary.Momentum.x;
                    if(is2D)
                        MomentumY[edge, n] = stateBoundary.Momentum.y;
                    if(is3D)
                        MomentumZ[edge, n] = stateBoundary.Momentum.z;
                    Energy[edge, n] = stateBoundary.Energy;
                }
            }



        }
    }
}
