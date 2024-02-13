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
    public class HLLCFlux_XDG_Interface : ILevelSetForm,  ISupportsJacobianComponent {

        #region ILevelSetForm members
        public int LevelSetIndex {
            get;
            private set;
        }

        public string PositiveSpecies {
            get;
            private set;
        }

        public string NegativeSpecies {
            get;
            private set;
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

        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            return this.bulkFlux.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB);
        }
            #endregion

        protected readonly HLLCFlux bulkFlux;

        public HLLCFlux_XDG_Interface(LevelSetTracker levelSetTracker, IBoundaryConditionMap boundaryConditionMap, Material material, FluxVariables flux, int component = int.MinValue, int levelSetIndex = 0, string posSpecies = "B", string negSpecies = "A") {
            //this.levelSetTracker = levelSetTracker;

            switch (flux) {
                case FluxVariables.Density:
                    this.bulkFlux = new HLLCDensityFlux(boundaryConditionMap, material);
                    break;
                case FluxVariables.Momentum:
                    this.bulkFlux = new HLLCMomentumFlux(boundaryConditionMap, component, material);
                    break;
                case FluxVariables.Energy:
                    this.bulkFlux = new HLLCEnergyFlux(boundaryConditionMap, material);
                    break;
                default:
                    throw new NotImplementedException();
            }

            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
        }
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if(SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");
            return new IEquationComponent[] {
                new LevelSetFormDifferentiator(this,SpatialDimension)
            };
        }
    }
}
