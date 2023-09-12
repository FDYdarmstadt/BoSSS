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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace XESF.Fluxes {
    class GodunovFlux_Wall : ILevelSetForm {

        #region ILevelSetForm members
        public int LevelSetIndex {
            get {
                return 0;
            }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }
        #endregion

        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            StateVector stateIn = new StateVector(material, density: 0.0, momentum: new Vector(0.0, 0.0), energy: 0.0);
            StateVector stateOut = new StateVector(uB, material);

            return this.bulkFlux.InnerEdgeFlux(null, int.MinValue, stateIn, stateOut, ref inp.Normal, int.MaxValue) * (0.0 - vB);
        }

        //private readonly LevelSetTracker levelSetTracker;

        private readonly GodunovFlux bulkFlux;

        private readonly Material material;

        public GodunovFlux_Wall(int SpatialDim, CompressibleControl control, IBoundaryConditionMap boundaryMap, Material material, FluxVariables flux, int component = int.MinValue) {
            //this.levelSetTracker = levelSetTracker;
            this.material = material;

            switch (flux) {
                case FluxVariables.Density:
                    this.bulkFlux = new GodunovFlux(control, boundaryMap, new EulerDensityComponent(), material);
                    break;
                case FluxVariables.Momentum:
                    this.bulkFlux = new GodunovFlux(control, boundaryMap, new EulerMomentumComponent(component, material.EquationOfState.HeatCapacityRatio, control.MachNumber, SpatialDim), material);
                    break;
                case FluxVariables.Energy:
                    this.bulkFlux = new GodunovFlux(control, boundaryMap, new EulerEnergyComponent(), material);
                    break;
                default:
                    throw new NotImplementedException();
            }
        }
    }
}
