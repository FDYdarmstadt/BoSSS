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
using log4net.Appender;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BUIDT.Fluxes
{
    /// <summary>
    /// S(pace) - T(ime) - Burgers - Upwind - Flux
    /// implements the weak form of the XDG discretization for the space time burgers equation using an smooth/unsmooth upwind flux
    /// </summary>
    class STBurgersUpwindFlux : ISpeciesFilter, IEdgeForm, IVolumeForm, ISupportsJacobianComponent
    {

        bool is_nf_smth;
        double s_alpha = 10;
        public Func<double[], double> DirichletBoundaryMap;

        public STBurgersUpwindFlux(string spc, bool is_nf_smth = false, double s_alpha = 10, Func<double[], double> dirichletBoundaryMap = null) {
            this.spc = spc;
            this.is_nf_smth = is_nf_smth;
            this.s_alpha = s_alpha;
            DirichletBoundaryMap = dirichletBoundaryMap;    
        }

        Vector FlowField(double[] x, double[] Uin, double[] Uout) {
            Vector u = new Vector(2);
            u.x = (Uin[0] + Uout[0]) * 2;
            u.y = 1;
            return u;
        }
        double SmoothedHeaviSide(double x) {
            return 1 / (1 + Math.Exp(-2 * s_alpha * x));
        }


        protected double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            Vector n = new Vector(2); n.x = normal[0]; n.y = normal[1];
            double[] Uout = new double[] { DirichletBoundaryMap(x) };
            var flowfield = FlowField(x, Uin, Uout);
            if(flowfield * n <= 0) { 
                Uout[0] = DirichletBoundaryMap(x);
            } else {
                Uout = Uin;
            }

            return InnerEdgeFlux(time, x, normal, Uin, Uout, jEdge);
        }

        protected double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge) {
            Vector n = new Vector(normal);
            var vel = FlowField(x, Uin, Uout);
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
            return ret;
        }

        protected void Flux(double time, double[] x, double[] U, double[] output) {
            output[0] = U[0] * U[0] * 0.5;
            output[1] = U[0];
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return InnerEdgeFlux(inp.time, inp.X, inp.Normal, _uIN, _uOUT, inp.iEdge) * (_vIN - _vOUT);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return BorderEdgeFlux(inp.time, inp.X, inp.Normal, inp.EdgeTag, _uA, inp.iEdge) * _vA;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double[] output = new double[2];
            Flux(cpv.time, cpv.Xglobal, U, output);
            return -1 * (output * (Vector)GradV);
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if(SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");

            return new IEquationComponent[] {
                new EdgeFormDifferentiator(this, SpatialDimension),
                new VolumeFormDifferentiator(this, SpatialDimension)
            };
        }
        public IList<string> ArgumentOrdering {
            get { return new string[] { "c" }; }
        }
        string spc = "";
        public string ValidSpecies => spc;

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public IList<string> ParameterOrdering => null;
    }

}