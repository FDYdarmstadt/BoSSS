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
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using ilPSP;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Symmetric Interior Penalty - discetization of the (positive) Laplace operator,
    /// i.e. \f$ + \text{div}( \nu \nabla u ) \f$.
    /// </summary>
    abstract public class SIPLaplace : IEdgeForm, IVolumeForm, IEquationComponentCoefficient, ISupportsJacobianComponent, IDGdegreeConstraint {

        /// <summary>
        /// no parameters in default implementation
        /// </summary>
        virtual public IList<string> ParameterOrdering { get { return null; } }

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="penalty_base">
        /// basic value for penalty parameter; in each cell it will be scaled 
        /// (in the order of) the inverse length scale and the squared polynomial degree
        /// </param>
        /// <param name="ArgumentVarName">
        /// the one and only string that is returned by the default implementation of <see cref="ArgumentOrdering"/>;
        /// </param>
        public SIPLaplace(double penalty_base, string ArgumentVarName) {
            m_penalty_base = penalty_base;
            m_ArgumentOrdering = new string[] { ArgumentVarName };
        }

        /// <summary>
        /// penalty parameter base multiplyer
        /// </summary>
        protected double m_penalty_base;

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


        /// <summary>
        /// update of penalty length scales.
        /// </summary>
        public virtual void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            
            double _D = cs.GrdDat.SpatialDimension;
            double _p = DomainDGdeg.Max();
            
            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty_deg = Math.Max(penalty_deg_tri, penalty_deg_sqr);
            
            this.LengthScales = cs.CellLengthScales;

            //this.LengthScales = ((Foundation.Grid.Classic.GridData)(cs.GrdDat)).Cells.CellLengthScale;
        }

        /// <summary>
        /// penalty degree multiplier
        /// </summary>
        protected double m_penalty_deg;
        
        /// <summary>
        /// Length scales used in <see cref="GetPenalty"/>
        /// </summary>
        protected MultidimensionalArray LengthScales;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        protected virtual double GetPenalty(int jCellIn, int jCellOut) {
            double cj_in = 1.0/LengthScales[jCellIn];
            
            double mu = m_penalty_base*m_penalty_deg * cj_in;
            if(jCellOut >= 0) {
                double cj_out = 1.0/LengthScales[jCellOut];
                mu = Math.Max(mu, m_penalty_base*m_penalty_deg * cj_out);
            }
            if(mu.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");

           
            return mu;
        }

        /// <summary>
        /// a little switch; turns everything off, except the penalty terms.
        /// </summary>
        protected double m_alpha = 1.0;


        /// <summary>
        /// Volume terms plus source terms for non-homogeneous boundary conditions
        /// </summary>
        virtual public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV);
            }
        }

        /// <summary>
        /// as induced by the SIP method
        /// </summary>
        virtual public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        /// <summary>
        /// as induced by the SIP method
        /// </summary>
        virtual public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }


        /// <summary>
        /// Volume integrand of the SIP
        /// </summary>
        public virtual double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int d = 0; d < cpv.D; d++)
                acc -= GradU[0, d] * GradV[d] * this.Nu(cpv.Xglobal, cpv.Parameters, cpv.jCell);
            acc *= this.m_alpha;
            return acc;
        }

        /// <summary>
        /// Integrand on interior mesh edges of the SIP
        /// </summary>
        virtual public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.GetPenalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);
            double nuB = this.Nu(inp.X, inp.Parameters_OUT, inp.jCellOut);


            for(int d = 0; d < inp.D; d++) {
                Acc += 0.5 * (nuA * _Grad_uA[0, d] + nuB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normal[d];  // consistency term
                Acc += 0.5 * (nuA * _Grad_vA[d] + nuB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normal[d];  // symmetry term
            }
            Acc *= this.m_alpha;

            double nuMax = (Math.Abs(nuA) > Math.Abs(nuB)) ? nuA : nuB;


            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * nuMax; // penalty term


            return Acc;

        }

        /// <summary>
        /// Integrand on boundary mesh edges of the SIP
        /// </summary>
        virtual public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if(this.IsDirichlet(ref inp)) {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);

                for(int d = 0; d < inp.D; d++) {
                    double nd = inp.Normal[d];
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

        /// <summary>
        /// Linear component - derivative is just this.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new[] { this };
        }

        /// <summary>
        /// checks that we have at least DG degree 1
        /// </summary>
        public bool IsValidDomainDegreeCombination(int[] DomainDegreesPerVariable, int CodomainDegree) {
            if (DomainDegreesPerVariable.Min() < 1)
                return false;
            if (CodomainDegree < 1)
                return false;
            
            return true;
        }
    }

}
