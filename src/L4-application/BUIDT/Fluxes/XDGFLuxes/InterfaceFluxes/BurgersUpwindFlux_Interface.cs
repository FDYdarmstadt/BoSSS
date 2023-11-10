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

namespace BUIDT.Fluxes {

    /// <summary>
    /// S(pace) - T(ime) - Burgers - Upwind - Flux - (over the) Interface
    /// implements the weak form of the XDG discretization for the space time burgers equation using an smooth/un-smooth upwind flux
    /// </summary>
    public class BurgersUpwindFlux_Interface : ILevelSetForm, ILevelSetEquationComponentCoefficient, ISupportsJacobianComponent {

        public double s_alpha;
        public bool is_nf_smth;

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

        Vector FlowField(double[] x, double[] Uin, double[] Uout) {
            Vector u = new Vector(2);
            u.x = (Uin[0] + Uout[0]) / 2;
            u.y = 1;
            return u;
        }
        double SmoothedHeaviSide(double x) {
            return 1 / (1 + Math.Exp(-2 * s_alpha * x));
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] Uin, double[] Uout, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            // Flux across the interface
            // Took Regular Flux
            Vector n = new Vector(2); n.x = inp.Normal.x; n.y = inp.Normal.y;
            var vel = FlowField(inp.X, Uin, Uout);
            double beta_n = vel * n;
            double ret;
            if(is_nf_smth) {
                ret = (0.5 * Uin[0] * Uin[0] * n.x + Uin[0] * n.y) * SmoothedHeaviSide(beta_n) + (0.5 * Uout[0] * Uout[0] * n.x + Uout[0] * n.y) * (1 - SmoothedHeaviSide(beta_n));
            } else {
                if(beta_n <= 0) {
                    ret = 0.5 * Uout[0] * Uout[0] * n.x + Uout[0] * n.y;
                } else {
                    ret = 0.5 * Uin[0] * Uin[0] * n.x + Uin[0] * n.y;
                }
            }
            ret *= (vA - vB);
            //Console.WriteLine($"Flux({Uin[0]}, {Uout[0]}) = {ret}");
            return ret;
        }



        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            return;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if(SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");
            return new IEquationComponent[] {
                new LevelSetFormDifferentiator(this,SpatialDimension)
            };
        }

        public bool IgnoreVectorizedImplementation => false;

        public IList<string> ParameterOrdering => null;




        //private readonly LevelSetTracker levelSetTracker;


        public BurgersUpwindFlux_Interface( bool is_nf_smth, double s_alpha, int levelSetIndex = 0, string posSpecies = "R", string negSpecies = "L") {
            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
            this.s_alpha = s_alpha;
            this.is_nf_smth = is_nf_smth;
        }
        public BurgersUpwindFlux_Interface() {

        }
    }
}
