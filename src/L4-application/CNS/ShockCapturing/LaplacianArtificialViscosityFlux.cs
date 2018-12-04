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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using ilPSP;
using System.Collections.Generic;

namespace CNS.ShockCapturing {

    class LaplacianArtificialViscosityFlux : SIPLaplace {

        private IBoundaryConditionMap boundaryMap;

        private ISpeciesMap speciesMap;

        //public LaplacianArtificialViscosityFlux(int order, MultidimensionalArray cj, string[] arguments, int affectedComponent, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap) :
        //    base((order + 1) * (order + CNSEnvironment.NumberOfDimensions) / (double)CNSEnvironment.NumberOfDimensions, cj, arguments, affectedComponent) {
        public LaplacianArtificialViscosityFlux(int order, MultidimensionalArray cj, Variable variable, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap) :
        base((order + 1) * (order + CNSEnvironment.NumberOfDimensions) / (double)CNSEnvironment.NumberOfDimensions, cj, variable) {
            this.boundaryMap = boundaryMap;
            this.speciesMap = speciesMap;
        }

        public override double Nu(double[] x, double[] parameter, int jCell) {
            return -1.0 * parameter[0];
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            return false;
        }

        //public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
        //    StateVector stateIn = new StateVector(_uA, speciesMap.GetMaterial(double.NaN));
        //    //StateVector boundaryState = boundaryMap.GetBoundaryState(inp.EdgeTag, inp.time, inp.X, inp.Normale, stateIn);
        //    //double[] _uB = boundaryState.ToArray();

        //    double[] _uB = stateIn.ToArray();

        //    double[,] _Grad_uB = _Grad_uA;

        //    double _vB = 0; // JA
        //    double[] _Grad_vB = _Grad_vA; // JA

        //    double[] paramsOut = inp.Parameters_IN; // Vielleicht AV innen = außen?
        //    //double[] paramsOut = new double[inp.Parameters_IN.Length]; // Vielleicht AV innen = außen?

        //    CommonParams inpInner = new CommonParams() {
        //        GridDat = inp.GridDat,
        //        iEdge = inp.iEdge,
        //        Normale = inp.Normale,
        //        Parameters_IN = inp.Parameters_IN,
        //        Parameters_OUT = paramsOut,
        //        time = inp.time,
        //        X = inp.X
        //    };

        //    return base.InnerEdgeForm(ref inpInner, _uA, _uB, _Grad_uA, _Grad_uB, _vA, _vB, _Grad_vA, _Grad_vB);
        //}

        //public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
        //    double Acc = 0.0;

        //    double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

        //    for (int d = 0; d < inp.D; d++)
        //        Acc += nuA * _Grad_uA[0, d] * _vA * inp.Normale[d]; // consistency term

        //    Acc *= this.m_alpha;

        //    return Acc;
        //}

        public override IList<string> ParameterOrdering {
            get {
                return new string[] { "artificialViscosity" };
            }
        }
    }
}
