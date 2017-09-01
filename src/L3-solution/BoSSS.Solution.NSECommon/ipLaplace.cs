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
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using ilPSP;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Interior Penalty - discetization of the Laplace operator
    /// </summary>
    abstract public class ipLaplace : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm {

        /// <summary>
        /// no parameters in default implementation
        /// </summary>
        virtual public IList<string> ParameterOrdering { get { return null; } }

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="penalty">
        /// basic value for penalty parameter: suggested to be <it>a*(p+1)*(p+D)/D</it>, where <it>p</it> is the polynomial degree of the involved DG field,
        /// <it>A</it> is the spatial dimension and <it>a</it> is in the order of 1, e.g. 0.9 or 1.2; see <see cref="GetPenalty"/>;
        /// </param>
        /// <param name="ArgumentVarName">
        /// the one and only string that is returned by the default implementation of <see cref="ArgumentOrdering"/>;
        /// </param>
        /// <param name="PenaltyLengthScales"></param>
        public ipLaplace(double penalty, MultidimensionalArray PenaltyLengthScales, string ArgumentVarName) {
            m_penalty = penalty;
            m_ArgumentOrdering = new string[] { ArgumentVarName };
            this.PenaltyLengthScales = PenaltyLengthScales;
        }

        /// <summary>
        /// penalty parameter (provided by ctor)
        /// </summary>
        private double m_penalty;

        /// <summary>
        /// Dirichlet boundary value
        /// </summary>
        virtual protected double g_Diri(ref Foundation.CommonParamsBnd inp) { return 0; }

        /// <summary>
        /// Neumann boundary value
        /// </summary>
        virtual protected double g_Neum(ref Foundation.CommonParamsBnd inp) { return 0; }


        string[] m_ArgumentOrdering;

        /// <summary>
        /// returns one argument variable, provided by the constructor
        /// </summary>
        virtual public IList<String> ArgumentOrdering { get { return m_ArgumentOrdering; } }

        /// <summary>
        /// diffusion coefficient, set to 1.0 per default;
        /// </summary>
        /// <param name="x">spatial position</param>
        /// <param name="p">parameter values</param>
        /// <param name="jCell">Cell index of the current cell</param>
        /// <returns></returns>
        virtual public double Nu(double[] x, double[] p, int jCell) {
            return 1.0;
        }

        /// <summary>
        /// true if the actual node (<paramref name="inp"/>, see <see cref="BoSSS.Foundation.InParams.X"/>),
        /// at which the flux is evaluated
        /// is in the Dirichlet domain
        /// </summary>
        /// <returns>
        /// if true, the actual node is assumed to be a Dirichlet boundary: <see cref="g_Diri"/>; <br/>
        /// if false, the actual node is assumed to be a Neumann boundary: <see cref="g_Neum"/>; <br/>
        /// </returns>
        protected abstract bool IsDirichlet(ref Foundation.CommonParamsBnd inp);

        protected  MultidimensionalArray PenaltyLengthScales;

       
        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        /// <param name="inp"></param>
        /// <returns></returns>
        protected virtual double GetPenalty(int jCellIn, int jCellOut) {
            double cj_in = PenaltyLengthScales[jCellIn];
            double mu = m_penalty* cj_in;
            if(jCellOut >= 0) {
                double cj_out = PenaltyLengthScales[jCellOut];
                mu = Math.Max(mu, m_penalty * cj_out);
            }

            return mu;
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

        virtual public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        public double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int d = 0; d < cpv.D; d++)
                acc -= GradU[0, d] * GradV[d] * this.Nu(cpv.Xglobal, cpv.Parameters, cpv.jCell) * this.m_alpha;
            return acc;
        }

        virtual public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.GetPenalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);
            double nuB = this.Nu(inp.X, inp.Parameters_OUT, inp.jCellOut);


            for(int d = 0; d < inp.D; d++) {
                Acc += 0.5 * (nuA * _Grad_uA[0, d] + nuB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                Acc += 0.5 * (nuA * _Grad_vA[d] + nuB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
            }
            Acc *= this.m_alpha;

            double nuMax = (Math.Abs(nuA) > Math.Abs(nuB)) ? nuA : nuB;


            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * nuMax; // penalty term
                        

            return Acc;

        }


        public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if(this.IsDirichlet(ref inp)) {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);

                for(int d = 0; d < inp.D; d++) {
                    double nd = inp.Normale[d];
                    Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                    Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                }
                Acc *= this.m_alpha;

                Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty

            } else {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }

    }

}
