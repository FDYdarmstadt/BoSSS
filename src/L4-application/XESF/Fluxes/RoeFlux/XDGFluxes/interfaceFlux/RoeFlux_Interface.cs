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
    public class RoeFlux_Interface : ILevelSetForm, IEquationComponentChecking, ISupportsJacobianComponent {

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
                return bulkFlux.ArgumentOrdering;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }
        #endregion
        public bool IgnoreVectorizedImplementation => false;

        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            return this.bulkFlux.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB);
        }

        private readonly LevelSetTracker levelSetTracker;

        protected readonly RoeBaseFlux bulkFlux;
        int D;

        public RoeFlux_Interface(LevelSetTracker levelSetTracker, IBoundaryConditionMap boundaryConditionMap, Material material, FluxVariables flux, int component, int levelSetIndex, string posSpecies, string negSpecies , double s_alpha, int D=2) {
            this.levelSetTracker = levelSetTracker;

            switch (flux) {
                case FluxVariables.Density:
                    this.bulkFlux = new RoeDensityFlux(boundaryConditionMap, material,s_alpha,D);
                    break;
                case FluxVariables.Momentum:
                    this.bulkFlux = new RoeMomentumFlux(boundaryConditionMap, component, material,s_alpha,D);
                    break;
                case FluxVariables.Energy:
                    this.bulkFlux = new RoeEnergyFlux(boundaryConditionMap, material, s_alpha,D);
                    break;
                default:
                    throw new NotImplementedException();
            }
            this.D = D;
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
