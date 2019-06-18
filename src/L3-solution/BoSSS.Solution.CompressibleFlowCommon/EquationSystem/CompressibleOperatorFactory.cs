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
using BoSSS.Foundation.Grid;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;

namespace BoSSS.Solution.CompressibleFlowCommon.EquationSystem {

    public static class CompressibleOperatorFactory {

        public static SpatialOperator BuildEulerOperator(IGridData gridData, CompressibleControl control) {

            // Boundary condition map
            Material material = control.GetMaterial();
            IBoundaryConditionMap boundaryMap = new BoundaryConditionMap(gridData, control, material);

            // Initialize operator
            SpatialOperator EulerOperator = new SpatialOperator(
                new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy },
                new string[] { },
                new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy },
                QuadOrderFunc.NonLinearWithoutParameters(2)
                );

            // Map fluxes
            EulerOperator.EquationComponents[CompressibleVariables.Density].Add(new OptimizedHLLCDensityFlux(boundaryMap, material));
            EulerOperator.EquationComponents[CompressibleVariables.Momentum.xComponent].Add(new OptimizedHLLCMomentumFlux(boundaryMap, 0, material));
            EulerOperator.EquationComponents[CompressibleVariables.Momentum.yComponent].Add(new OptimizedHLLCMomentumFlux(boundaryMap, 1, material));
            EulerOperator.EquationComponents[CompressibleVariables.Energy].Add(new OptimizedHLLCEnergyFlux(boundaryMap, material));

            EulerOperator.Commit();

            return EulerOperator;
        }
    }
}
