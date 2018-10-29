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

using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using CNS.MaterialProperty;
using ilPSP;

namespace CNS.Convection {

    /// <summary>
    /// Base class for optimized versions of the HLLC flux
    /// </summary>
    public abstract class OptimizedHLLCFlux : INonlinearFlux {

        /// <summary>
        /// <see cref="OptimizedHLLCDensityFlux.OptimizedHLLCDensityFlux"/>
        /// </summary>
        protected readonly CNSControl config;

        /// <summary>
        /// <see cref="OptimizedHLLCDensityFlux.OptimizedHLLCDensityFlux"/>
        /// </summary>
        protected readonly ISpeciesMap speciesMap;

        /// <summary>
        /// <see cref="OptimizedHLLCDensityFlux.OptimizedHLLCDensityFlux"/>
        /// </summary>
        protected readonly IBoundaryConditionMap boundaryMap;

        /// <summary>
        /// Constructs a new flux
        /// </summary>
        /// <param name="config">
        /// Configuration options
        /// </param>
        /// <param name="speciesMap">
        /// Species map. Only support ideal gas in the entire domain.
        /// </param>
        /// <param name="boundaryMap">
        /// Mapping for boundary conditions
        /// </param>
        public OptimizedHLLCFlux(CNSControl config, ISpeciesMap speciesMap, IBoundaryConditionMap boundaryMap) {
            this.config = config;
            this.speciesMap = speciesMap;
            this.boundaryMap = boundaryMap;
        }

        #region INonlinearFlux Members

        /// <summary>
        /// <see cref="INonlinearFlux.BorderEdgeFlux"/>
        /// </summary>
        /// <param name="time"></param>
        /// <param name="jEdge"></param>
        /// <param name="x"></param>
        /// <param name="normal"></param>
        /// <param name="normalFlipped"></param>
        /// <param name="EdgeTags"></param>
        /// <param name="EdgeTagsOffset"></param>
        /// <param name="Uin"></param>
        /// <param name="Offset"></param>
        /// <param name="Lenght"></param>
        /// <param name="Output"></param>
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
            int D = CNSEnvironment.NumberOfDimensions;
            double sign = normalFlipped ? -1.0 : 1.0;

            MultidimensionalArray[] Uout = new MultidimensionalArray[Uin.Length];
            for (int i = 0; i < Uin.Length; i++) {
                Uout[i] = MultidimensionalArray.Create(Uin[i].GetLength(0), Uin[i].GetLength(1));
            }

            BoSSS.Solution.CompressibleFlowCommon.MaterialProperty.Material material = speciesMap.GetMaterial(double.NaN);
            for (int e = 0; e < Lenght; e++) {
                int edge = e + Offset;
                for (int n = 0; n < NoOfNodes; n++) {
                    double[] xLocal = new double[D];
                    double[] normalLocal = new double[D];
                    for (int d = 0; d < D; d++) {
                        xLocal[d] = x[edge, n, d];
                        normalLocal[d] = normal[edge, n, d] * sign;
                    }

                    StateVector stateIn = new StateVector(material, Uin, edge, n, D);
                    StateVector stateBoundary = boundaryMap.GetBoundaryState(
                        EdgeTags[e + EdgeTagsOffset], time, xLocal, normalLocal, stateIn);

                    Uout[0][edge, n] = stateBoundary.Density;
                    for (int d = 0; d < D; d++) {
                        Uout[d + 1][edge, n] = stateBoundary.Momentum[d];
                    }
                    Uout[D + 1][edge, n] = stateBoundary.Energy;
                }
            }

            InnerEdgeFlux(time, jEdge, x, normal, Uin, Uout, Offset, Lenght, Output);
        }

        /// <summary>
        /// <see cref="INonlinearFlux.InnerEdgeFlux"/>
        /// </summary>
        /// <param name="time"></param>
        /// <param name="jEdge"></param>
        /// <param name="x"></param>
        /// <param name="normal"></param>
        /// <param name="Uin"></param>
        /// <param name="Uout"></param>
        /// <param name="Offset"></param>
        /// <param name="Lenght"></param>
        /// <param name="Output"></param>
        public abstract void InnerEdgeFlux(double time, int jEdge, MultidimensionalArray x, MultidimensionalArray normal, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, int Offset, int Lenght, MultidimensionalArray Output);

        /// <summary>
        /// <see cref="INonlinearFlux.Flux"/>
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="U"></param>
        /// <param name="Offset"></param>
        /// <param name="Length"></param>
        /// <param name="Output"></param>
        public abstract void Flux(double time, MultidimensionalArray x, ilPSP.MultidimensionalArray[] U, int Offset, int Length, MultidimensionalArray Output);

        #endregion

        #region IEquationComponent Members

        /// <summary>
        /// <see cref="CNSEnvironment.PrimalArgumentOrdering"/>
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return CNSEnvironment.PrimalArgumentOrdering;
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
