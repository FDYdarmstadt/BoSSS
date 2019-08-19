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
using ilPSP;
using System.Collections.Generic;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockCapturing {

    public class LaplacianArtificialViscosityFlux : SIPLaplace {

        private readonly IBoundaryConditionMap boundaryMap;

        private readonly ISpeciesMap speciesMap;

        public LaplacianArtificialViscosityFlux(int order, MultidimensionalArray cj, Variable variable, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap) :
        base((order + 1) * (order + CompressibleEnvironment.NumberOfDimensions) / (double)CompressibleEnvironment.NumberOfDimensions, cj, variable) {
            this.boundaryMap = boundaryMap;
            this.speciesMap = speciesMap;
        }

        public override double Nu(double[] x, double[] parameter, int jCell) {
            return -1.0 * parameter[0];
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            return false;
        }

        public override IList<string> ParameterOrdering {
            get {
                return new string[] { "artificialViscosity" };
            }
        }
    }
}
