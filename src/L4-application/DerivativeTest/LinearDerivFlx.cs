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

namespace BoSSS.Application.DerivativeTest {
    
    
    public class LinearDerivFlx : LinearFlux {
        public LinearDerivFlx(int __d) {
            d = __d;
        }

        int d;

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            //return Uin[0]*inp.Normale[d];
            return inp.Parameters_IN[0] * inp.Normale[d];
        }

        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            return 0.5*(Uin[0] + Uout[0])*inp.Normale[d];
        }

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            Debug.Assert(output.Length == inp.Xglobal.Length);
            Array.Clear(output, 0, output.Length);
            output[d] = U[0];
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { "u" }; 
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return new string[] { "uBnd" }; 
            }
        }
    }

    /// <summary>
    /// Interior Penalty - discetization of the Laplace operator
    /// </summary>
    class ipLaplace : BoSSS.Foundation.IVolumeForm, BoSSS.Foundation.IEdgeForm {

        /// <summary>
        /// no parameters in default implementation
        /// </summary>
        virtual public IList<string> ParameterOrdering { get { return null; } }

        /// <summary>
        /// ctor
        /// </summary>
        public ipLaplace() {
            m_ArgumentOrdering = new string[] { "u" };
        }
        
        string[] m_ArgumentOrdering;

        /// <summary>
        /// returns one argument variable, provided by the constructor
        /// </summary>
        public IList<String> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
            }
        }

        /// <summary>
        /// diffusion coefficient, set to 1.0 per default;
        /// </summary>
        /// <param name="x">spatial position</param>
        /// <param name="p">parameter values</param>
        /// <returns></returns>
        double Nu(double[] x, double[] p) {
            return 1.0;
        }



        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        /// <param name="inp"></param>
        /// <returns></returns>
        protected virtual double mu(int jCellIn, int jCellOut) {
            /*
            double cj_in = cj[jCellIn];
            double mu = m_penalty * cj_in;
            if(jCellOut >= 0) {
                double cj_out = cj[jCellOut];
                mu = Math.Max(mu, m_penalty * cj_out);
            }
            //if(!rem) {
            //    rem = true;
            //    Console.WriteLine("REM: fake penalty parameter");
            //}
            //mu = 30;
            return mu;
            */
            return 0.0;
        }


        /// <summary>
        /// a little switch...
        /// </summary>
        protected double m_alpha = 1.0;


        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV);
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        public double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int d = 0; d < cpv.D; d++)
                acc -= GradU[0, d] * GradV[d] * this.Nu(cpv.Xglobal, cpv.Parameters) * this.m_alpha;
            return acc;
        }

        virtual public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.mu(inp.jCellIn, inp.jCellOut);
            double muA = this.Nu(inp.X, inp.Parameters_IN);
            double muB = this.Nu(inp.X, inp.Parameters_OUT);

            for(int d = 0; d < inp.D; d++) {
                Acc += 0.5 * (muA * _Grad_uA[0, d] + muB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                Acc += 0.5 * (muA * _Grad_vA[d] + muB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
            }
            Acc *= this.m_alpha;

            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;


            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * muMax; // penalty term

            return Acc;
        }

        public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.mu(inp.jCellIn, -1);
            double muA = this.Nu(inp.X, inp.Parameters_IN);

            {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = _uA[0];

                for(int d = 0; d < inp.D; d++) {
                    double nd = inp.Normale[d];
                    Acc += (muA * _Grad_uA[0, d]) * (_vA) * nd;
                    Acc += (muA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                }
                Acc *= this.m_alpha;

                Acc -= muA * (_uA[0] - g_D) * (_vA - 0) * pnlty;

            } 
            return Acc;
        }

    }

}
