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

namespace SAIDT.Fluxes {
    class ScalarAdvectionUpwindFlux : ISpeciesFilter, IEdgeForm, IVolumeForm {
        public double cL;
        public double cR;
        public double x0;

        public Func<double[], double> DirichletBoundaryMap;
        public Func<double, double> FlowFunc;

        public ScalarAdvectionUpwindFlux(string spc, double x0, double cL, double cR, Func<double, double> FlowFunc, Func<double[], double> bndrymap) {
            this.cL = cL;
            this.cR = cR;
            this.spc = spc;
            this.x0 = x0;
            this.FlowFunc = FlowFunc;
            this.DirichletBoundaryMap = bndrymap;
        }
        public Vector FlowField(double[] X) {
            var flowfield = new Vector(2);
            flowfield.x = FlowFunc(X[1]);
            flowfield.y = 1;
            return flowfield;
        }

        protected double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            Vector n = new Vector(2); n.x = normal[0]; n.y = normal[1];
            var flowfield = FlowField(x);
            double[] Uout = new double[1];
            if(flowfield * n < 0) { // inflow
                // if(x[0] > x0) {
                //     Uout[0] = cR;
                // } else {
                //     Uout[0] = cL;
                // }
                Uout[0] = DirichletBoundaryMap(x);
            } else {
                Uout = Uin;
            }
            return InnerEdgeFlux(time, x, normal, Uin, Uout, jEdge);
            // }else{
            //     if (x[0]<0.25){
            //         return vel*n;
            //     }else{
            //         return 0;
            //     }
            // }

            // var vel = FlowField(x, Uin, Uout);
            // return (vel * Uout[0]) * n;
        }

        protected double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge) {
            Vector n = new Vector(normal);
            var flowfield = FlowField(x);
            if(flowfield * n > 0)
                return (flowfield * Uin[0]) * n;
            else
                return (flowfield * Uout[0]) * n;
        }

        protected void Flux(double time, double[] x, double[] U, double[] output) {
            output[0] = FlowFunc(x[1]) * U[0];
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

        /// <summary>
        /// 
        /// </summary>
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