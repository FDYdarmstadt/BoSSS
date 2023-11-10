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

namespace SAIDT.Fluxes
{
    public class ScalarAdvectionUpwindFlux_Interface : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        private LevelSetTracker lsTrk;
        private double a;
        public Func<double, double> FlowFunc;

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
                return new string[] { "c" };
            }
        }

        public Vector FlowField(double[] X) {
            var flowfield = new Vector(2);
            flowfield.x = FlowFunc(X[1]);
            flowfield.y = 1;
            return flowfield;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            // Flux across the interface
            // Took Regular Flux
            Vector n = new Vector(2); n.x = inp.Normal.x; n.y = inp.Normal.y;
            var vel = FlowField(inp.X);

            if(vel * n > 0)
                return (vel * uA[0]) * n * (vA - vB);
            else
                return (vel * uB[0]) * n * (vA - vB);

            //throw new NotSupportedException("Has to be checked again.");
        }



        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            return;
        }


        public bool IgnoreVectorizedImplementation => false;

        public IList<string> ParameterOrdering => null;




        //private readonly LevelSetTracker levelSetTracker;

        public ScalarAdvectionUpwindFlux_Interface(LevelSetTracker levelSetTracker, int levelSetIndex = 0, string posSpecies = "R", string negSpecies = "L") {
            //this.levelSetTracker = levelSetTracker;

            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
        }

        public ScalarAdvectionUpwindFlux_Interface(LevelSetTracker lsTrk, Func<double, double> FlowFunc, int levelSetIndex = 0, string posSpecies = "R", string negSpecies = "L") {
            this.lsTrk = lsTrk;
            this.FlowFunc = FlowFunc;
            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
        }

        public ScalarAdvectionUpwindFlux_Interface(double a) {
            this.a = a;
        }
    }

}
