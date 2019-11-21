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


        public double LevelSetForm(ref CommonParamsLs inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            // Convert CommonParamsLs to CommonParams
            CommonParams commonParams;
            commonParams.Normal = inp.Normal;
            commonParams.X = inp.X;
            commonParams.Parameters_IN = inp.ParamsNeg;
            commonParams.Parameters_OUT = inp.ParamsPos;
            commonParams.iEdge = int.MaxValue; // Alternative: use cell index as edge index --> Might this cause problems?
            commonParams.GridDat = this.levelSetTracker.GridDat;
            commonParams.time = inp.time;

            return this.bulkFlux.InnerEdgeFormHack(ref commonParams, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB);
        }

        public void LevelSetForm_GradV(LevSetIntParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_GradV) {
            // Convert LevSetIntParams to EdgeFormParams
            EdgeFormParams edgeParams;
            edgeParams.e0 = inp.i0;
            edgeParams.Len = inp.Len;
            edgeParams.GridDat = null;
            edgeParams.time = inp.time;
            edgeParams.ParameterVars_IN = inp.ParamsNeg;
            edgeParams.ParameterVars_OUT = inp.ParamsPos;
            edgeParams.Normals = inp.Normal;
            edgeParams.NodesGlobal = inp.X;

            // Split coeff array into species A and B
            MultidimensionalArray fin = Koeff_GradV.ExtractSubArrayShallow(-1, -1, 0, -1);
            MultidimensionalArray fout = Koeff_GradV.ExtractSubArrayShallow(-1, -1, 1, -1);

            this.bulkFlux.InternalEdge_GradV(ref edgeParams, uA, uB, Grad_uA, Grad_uB, fin, fout);
        }

        public void LevelSetForm_V(LevSetIntParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_V) {
            // Convert LevSetIntParams to EdgeFormParams
            EdgeFormParams edgeParams;
            edgeParams.e0 = inp.i0;
            edgeParams.Len = inp.Len;
            edgeParams.GridDat = null;
            edgeParams.time = inp.time;
            edgeParams.ParameterVars_IN = inp.ParamsNeg;
            edgeParams.ParameterVars_OUT = inp.ParamsPos;
            edgeParams.Normals = inp.Normal;
            edgeParams.NodesGlobal = inp.X;

            // Split coeff array into species A and B
            MultidimensionalArray fin = Koeff_V.ExtractSubArrayShallow(-1, -1, 0);
            MultidimensionalArray fout = Koeff_V.ExtractSubArrayShallow(-1, -1, 1);

            this.bulkFlux.InternalEdge_V(ref edgeParams, uA, uB, Grad_uA, Grad_uB, fin, fout);
        }
    }
}
