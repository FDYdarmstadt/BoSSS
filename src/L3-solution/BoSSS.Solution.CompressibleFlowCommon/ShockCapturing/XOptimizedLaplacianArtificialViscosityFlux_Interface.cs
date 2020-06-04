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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockCapturing {

    public class XOptimizedLaplacianArtificialViscosityFlux_Interface : INonlinLevelSetForm_V, INonlinLevelSetForm_GradV {

        private readonly LevelSetTracker levelSetTracker;

        private readonly XOptimizedLaplacianArtificialViscosityFlux bulkFlux;

        public XOptimizedLaplacianArtificialViscosityFlux_Interface(LevelSetTracker levelSetTracker, string ArgumentVarName, double penaltySafetyFactor, double penaltyFactor, Dictionary<SpeciesId, MultidimensionalArray> inverseLengthScales) {
            this.levelSetTracker = levelSetTracker;
            this.ArgumentOrdering = new string[] { ArgumentVarName };
            this.bulkFlux = new XOptimizedLaplacianArtificialViscosityFlux(null, levelSetTracker, ArgumentVarName, penaltySafetyFactor, penaltyFactor, inverseLengthScales);
        }

        public int LevelSetIndex {
            get {
                return 0;
            }
        }

        public SpeciesId PositiveSpecies {
            get {
                return this.levelSetTracker.GetSpeciesId("B");
            }
        }

        public SpeciesId NegativeSpecies {
            get {
                return this.levelSetTracker.GetSpeciesId("A");
            }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ArgumentOrdering {
            get;
            private set;
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { "artificialViscosity" };
            }
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            return this.bulkFlux.InnerEdgeFormHack(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB);
        }

        public void NonlinInternalEdge_GradV(ref EdgeFormParams edgeParams, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_GradVin, MultidimensionalArray Koeff_GradVot) {
            this.bulkFlux.InternalEdge_GradV(ref edgeParams, uA, uB, Grad_uA, Grad_uB, Koeff_GradVin, Koeff_GradVot);
        }

        public void NonlinInternalEdge_V(ref EdgeFormParams edgeParams, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_Vin, MultidimensionalArray Koeff_Vot) {
            this.bulkFlux.InternalEdge_V(ref edgeParams, uA, uB, Grad_uA, Grad_uB, Koeff_Vin, Koeff_Vot);
        }
    }
}
