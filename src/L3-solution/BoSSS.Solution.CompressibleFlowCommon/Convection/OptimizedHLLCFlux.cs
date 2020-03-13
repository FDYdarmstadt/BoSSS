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
using System.Diagnostics;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;

namespace BoSSS.Solution.CompressibleFlowCommon.Convection {

    /// <summary>
    /// Base class for optimized versions of the HLLC flux
    /// </summary>
    public abstract class OptimizedHLLCFlux : INonlinearFlux, IEquationComponentSpeciesNotification {

        /// <summary>
        /// <see cref="OptimizedHLLCDensityFlux.OptimizedHLLCDensityFlux"/>
        /// </summary>
        protected readonly IBoundaryConditionMap boundaryMap;

        protected readonly Material material;

        private string speciesName;

        /// <summary>
        /// Constructs a new flux
        /// </summary>
        /// <param name="boundaryMap">
        /// Mapping for boundary conditions
        /// </param>
        public OptimizedHLLCFlux(IBoundaryConditionMap boundaryMap, Material material) {
            this.boundaryMap = boundaryMap;
            this.material = material;
        }

        #region INonlinearFlux Members

        /// <summary>
        /// <see cref="INonlinearFlux.BorderEdgeFlux"/>
        /// </summary>
        public void BorderEdgeFlux(
            double time,
            int jEdge,
            MultidimensionalArray x,
            MultidimensionalArray normal,
            bool normalFlipped,
            byte[] EdgeTags,
            int EdgeTagsOffset,
            MultidimensionalArray[] Uin,
            int Offset,
            int Lenght,
            MultidimensionalArray Output) {
            int NoOfNodes = Uin[0].GetLength(1);
            int D = CompressibleEnvironment.NumberOfDimensions;

            MultidimensionalArray[] Uout = new MultidimensionalArray[Uin.Length];
            for (int i = 0; i < Uin.Length; i++) {
                Uout[i] = MultidimensionalArray.Create(Uin[i].GetLength(0), Uin[i].GetLength(1));
            }

            var xdgBoudaryMap = this.boundaryMap as XDGCompressibleBoundaryCondMap;
            var boundaryMap = this.boundaryMap as CompressibleBoundaryCondMap;
            if (xdgBoudaryMap == null && boundaryMap == null)
                throw new NotSupportedException("This type of boundary condition map is not supported.");
            Vector xLocal = new Vector(D);
            Vector normalLocal = new Vector(D);

            bool isWall = false;

            // Loop over edges
            for (int e = 0; e < Lenght; e++) {
                byte EdgeTag = EdgeTags[e + EdgeTagsOffset];

                // Sweep until the boundary condition changes
                int __L = 1;
                for (; e + __L < Lenght; __L++) {
                    byte _EdgeTag = EdgeTags[e + __L + EdgeTagsOffset];
                    if (EdgeTag != _EdgeTag)
                        break;
                }

                // Get boundary condition on this edge
                BoundaryCondition boundaryCondition;
                if (xdgBoudaryMap != null)
                    boundaryCondition = xdgBoudaryMap.GetBoundaryConditionForSpecies(EdgeTags[e + EdgeTagsOffset], this.speciesName);
                else
                    boundaryCondition = boundaryMap.GetBoundaryCondition(EdgeTags[e + EdgeTagsOffset]);

                // Call vectorized state evaluation 
                int edge = e + Offset;
                boundaryCondition.GetBoundaryState(Uout, time, x, normal, Uin, edge, __L, normalFlipped, material);
                e += __L - 1;

                if (boundaryCondition is AdiabaticSlipWall) {
                    isWall = true;
                }
            }

            if (!isWall) {
                InnerEdgeFlux(time, jEdge, x, normal, Uin, Uout, Offset, Lenght, Output);
            }

            //if (writeFluxes) {
            //    Console.WriteLine(String.Format("Flux at ({0:0.00000000}, {1:0.00000000}) = {2:0.00000000} \t normal = ({3:0.00000000}, {4:0.00000000})", x[0, 0, 0], x[0, 0, 1], Output[0, 0], normal[0, 0, 0], normal[0, 0, 1]));
            //    //Console.WriteLine("hallo");
            //}
        }

        /// <summary>
        /// <see cref="INonlinearFlux.InnerEdgeFlux"/>
        /// </summary>
        public abstract void InnerEdgeFlux(double time, int jEdge, MultidimensionalArray x, MultidimensionalArray normal, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, int Offset, int Length, MultidimensionalArray Output);

        /// <summary>
        /// <see cref="INonlinearFlux.Flux"/>
        /// </summary>
        public abstract void Flux(double time, MultidimensionalArray x, ilPSP.MultidimensionalArray[] U, int Offset, int Length, MultidimensionalArray Output);

        public void SetParameter(string speciesName, SpeciesId SpcId) {
            this.speciesName = speciesName;
        }

        #endregion

        #region IEquationComponent Members

        /// <summary>
        /// <see cref="CompressibleEnvironment.PrimalArgumentOrdering"/>
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return CompressibleEnvironment.PrimalArgumentOrdering;
            }
        }

        /// <summary>
        /// Empty (i.e., no parameters are used)
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        #endregion
    }
}
