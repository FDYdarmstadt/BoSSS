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
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockCapturing {

    public class LaplacianArtificialViscosityFlux : SIPLaplace {

        private readonly BoundaryCondMap<XDGHeatBcType> boundaryCondMap;

        /// <summary>
        /// Implements the positive Laplace operator, inherits from <see cref="SIPLaplace"/>
        /// </summary>
        /// <param name="boundaryCondMap">Information about boundary conditions</param>
        /// <param name="penaltySafteyFactor">A user definded factor, typically in the range of 3.0 to 5.0</param>
        /// <param name="penaltyFactor">A factor based on the grid type (tetras, quads, etc.)</param>
        /// <param name="lengthScales">A cell length scale</param>
        /// <param name="argumentName">The variable where the operator acts on</param>
        public LaplacianArtificialViscosityFlux(BoundaryCondMap<XDGHeatBcType> boundaryCondMap, double penaltySafteyFactor, double penaltyFactor, MultidimensionalArray lengthScales, string argumentName) :
              base(penaltySafteyFactor * penaltyFactor, lengthScales, argumentName) {
            this.boundaryCondMap = boundaryCondMap;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            throw new NotSupportedException("I had to implement this...");
        }

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            XDGHeatBcType edgeType = this.boundaryCondMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case XDGHeatBcType.Dirichlet:
                    Func<double[], double, double> dirichletFunction = this.boundaryCondMap.bndFunction["u"][inp.EdgeTag];
                    double g_D = dirichletFunction(inp.X, inp.time);

                    for (int d = 0; d < inp.D; d++) {
                        double nd = inp.Normale[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                        Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                    }
                    Acc *= this.m_alpha;

                    Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty
                    break;

                case XDGHeatBcType.ZeroNeumann:
                    double g_N = 0.0;

                    Acc += nuA * g_N * _vA * this.m_alpha;
                    break;

                default:
                    throw new NotSupportedException();
            }

            return Acc;
        }
    }
}
