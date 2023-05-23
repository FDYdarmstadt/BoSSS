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

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Solution.Utils;
using BoSSS.Foundation;
using System.Diagnostics;
using ilPSP;
using BoSSS.Foundation.XDG;

namespace BoSSS.Application.DerivativeTest_XDG_Enriched {
    
    class UpwindFlux_XDG_Interface : ILevelSetForm, ILevelSetEquationComponentCoefficient
    {

        private LevelSetTracker lsTrk;
        private double a;
        public Func<double, double> FlowFunc;

        public int LevelSetIndex
        {
            get;
            private set;
        }

        public string PositiveSpecies
        {
            get;
            private set;
        }

        public string NegativeSpecies
        {
            get;
            private set;
        }

        public TermActivationFlags LevelSetTerms
        {
            get
            {
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ArgumentOrdering
        {
            get
            {
                return new string[] { "c" };
            }
        }

        public Vector FlowField(double[] X)
        {
            var flowfield = new Vector(2);
            flowfield.x = FlowFunc(X[1]);
            flowfield.y = 1;
            return flowfield;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB)
        {
            // Flux across the interface
            // Took Regular Flux
            Vector n = new Vector(2); n.x = inp.Normal.x; n.y = inp.Normal.y;
            var vel = FlowField(inp.X);

            if (vel * n > 0)
                return (vel * uA[0]) * n * (vA - vB);
            else
                return (vel * uB[0]) * n * (vA - vB);

            //throw new NotSupportedException("Has to be checked again.");
        }



        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg)
        {
            return;
        }


        public bool IgnoreVectorizedImplementation => false;

        public IList<string> ParameterOrdering => null;




        //private readonly LevelSetTracker levelSetTracker;

        public UpwindFlux_XDG_Interface(LevelSetTracker levelSetTracker, int levelSetIndex = 0, string posSpecies = "B", string negSpecies = "A")
        {
            //this.levelSetTracker = levelSetTracker;

            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
        }

        public UpwindFlux_XDG_Interface(LevelSetTracker lsTrk, Func<double, double> FlowFunc, int levelSetIndex = 0, string posSpecies = "B", string negSpecies = "A")
        {
            this.lsTrk = lsTrk;
            this.FlowFunc = FlowFunc;
            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
        }

        public UpwindFlux_XDG_Interface(double a)
        {
            this.a = a;
        }
    }

    class ScalarTransportFlux : ISpeciesFilter, IEdgeForm, IVolumeForm
    {
        public double cL;
        public double cR;
        public double x0;

        public Func<double, double> FlowFunc;

        public ScalarTransportFlux(string spc, double x0, double cL, double cR, Func<double, double> FlowFunc)
        {
            this.cL = cL;
            this.cR = cR;
            this.spc = spc;
            this.x0 = x0;
            this.FlowFunc = FlowFunc;
        }
        public Vector FlowField(double[] X)
        {
            var flowfield = new Vector(2);
            flowfield.x = FlowFunc(X[1]);
            flowfield.y = 1;
            return flowfield;
        }

        protected double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge)
        {
            Vector n = new Vector(2); n.x = normal[0]; n.y = normal[1];
            var flowfield = FlowField(x);
            double[] Uout = new double[1];
            if (flowfield * n < 0)
            { // inflow
                if (x[0] > x0)
                {
                    Uout[0] = cR;
                }
                else
                {
                    Uout[0] = cL;
                }
            }
            else
            {
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

        protected double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge)
        {
            Vector n = new Vector(normal);
            var flowfield = FlowField(x);
            if (flowfield * n > 0)
                return (flowfield * Uin[0]) * n;

            else
                return (flowfield * Uout[0]) * n;
        }

        protected void Flux(double time, double[] x, double[] U, double[] output)
        {
            var flowfield = FlowField(x);
            Vector o;
            o = flowfield * U[0];
            output[0] = o.x;
            output[1] = o.y;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT)
        {
            return InnerEdgeFlux(inp.time, inp.X, inp.Normal, _uIN, _uOUT, inp.iEdge) * (_vIN - _vOUT);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
        {
            return BorderEdgeFlux(inp.time, inp.X, inp.Normal, inp.EdgeTag, _uA, inp.iEdge) * _vA;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double[] output = new double[2];
            Flux(cpv.time, cpv.Xglobal, U, output);
            return -1 * (output * (Vector)GradV);
        }

        /// <summary>
        /// 
        /// </summary>
        public IList<string> ArgumentOrdering
        {
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
