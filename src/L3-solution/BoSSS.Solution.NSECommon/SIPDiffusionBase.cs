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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation;
using ilPSP;

namespace BoSSS.Solution.NSECommon
{

    /// <summary>
    /// SIP discretization of diffusion operators for scalar transport equations (i.e. species mass transport and temperature). Analogous to swipViscosity_Term1.
    /// </summary>
    public abstract class SIPDiffusionBase : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm, IEquationComponentCoefficient {
 
        /// <summary>
        /// The Function in \nabla \dot (Diffusivity \nabla u), e.g. heat conductivity or diffusion coefficient
        /// </summary>
        protected abstract double Diffusivity(double[] Parameters, double[,] GradU , Vector NodeCoordinates);

        protected double PenaltyBase;
        protected Func<double[], double, double>[] ArgumentFunction;

        /// <summary>
        /// Ctor for Species transport equations
        /// </summary>
        /// <param name="Coefficient">Coefficient Function in \nabla \dot (Coefficient \nabla u)</param>
        /// <param name="PenaltyBase">C.f. Calculation of SIP penalty base, cf. Chapter 3 in 
        /// K. Hillewaert, “Development of the discontinuous Galerkin method for high-resolution, large scale CFD and acoustics in industrial geometries”,
        /// Université catholique de Louvain, 2013.</param>
        /// <param name="BcMap">Boundary condition map</param>
        /// <param name="Argument">The argument of the flux. Must be compatible with the DiffusionMode.</param>
        /// <param name="PenaltyLengthScales"></param>
        protected SIPDiffusionBase(double PenaltyBase, MultidimensionalArray PenaltyLengthScales, bool ParametersOK = false, int speciesIndex = 0) {
            this.PenaltyBase = PenaltyBase;
            this.cj = PenaltyLengthScales;
            this.prmsOK = ParametersOK;
            this.i = speciesIndex;
        }

        /// <summary>
        /// Ctor for Species transport equations
        /// </summary>
        /// <param name="Coefficient">Coefficient Function in \nabla \dot (Coefficient \nabla u)</param>
        /// <param name="PenaltyBase">C.f. Calculation of SIP penalty base, cf. Chapter 3 in 
        /// K. Hillewaert, “Development of the discontinuous Galerkin method for high-resolution, large scale CFD and acoustics in industrial geometries”,
        /// Université catholique de Louvain, 2013.</param>
        /// <param name="BcMap">Boundary condition map</param>
        /// <param name="Argument">The argument of the flux. Must be compatible with the DiffusionMode.</param>
        /// <param name="PenaltyLengthScales"></param>
        protected SIPDiffusionBase(double PenaltyBase, bool ParametersOK = false, int speciesIndex = 0) {
            this.PenaltyBase = PenaltyBase;
            this.prmsOK = ParametersOK;
            this.i = speciesIndex;
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV);
                //return TermActivationFlags.AllOn;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
                //return TermActivationFlags.AllOn;
            }
        }

        public TermActivationFlags VolTerms {
            get {
                //return TermActivationFlags.AllOn;
                return TermActivationFlags.GradUxGradV | TermActivationFlags.UxGradV; // should be changed. The term UxGradV is not allways needed
            }
        }


        /// <summary>
        /// really really ugly hack for defining if the <see cref="Diffusivity(double[])"/> method should take arguments or variables as argument
        /// </summary>
        bool prmsOK;

        /// <summary>
        /// Index that indicates which element of the argumentOrdering is being used for calculating the forms. Default is 0
        /// </summary>
        protected int i;

        protected MultidimensionalArray cj;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        /// <returns></returns>
        protected double GetPenalty(int jCellIn, int jCellOut) {

            double penaltySizeFactor = 1;
            if (cj != null) {
                double penaltySizeFactor_A = 1.0 / cj[jCellIn];
                double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0;
                Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
                Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
                Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
                Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
                penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);
            }

            Debug.Assert(!double.IsInfinity(m_penalty));
            Debug.Assert(!double.IsInfinity(m_penalty));

            double µ = penaltySizeFactor * m_penalty * PenaltyBase;
            if(µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        /// <summary>
        /// spatial dimension
        /// </summary>
        protected int m_D;

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        protected double m_penalty;

        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        /// <param name="cs"></param>
        /// <param name="DomainDGdeg"></param>
        /// <param name="TestDGdeg"></param>
        public virtual void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            m_D = cs.GrdDat.SpatialDimension;
            double _D = m_D;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = cs.CellLengthScales;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;
            double pnlty = GetPenalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);

            double DiffusivityA;
            double DiffusivityB;
            double DiffusivityMax;
            double[] difusivityArguments_IN = prmsOK ? inp.Parameters_IN : _uA;
            double[] difusivityArguments_OUT = prmsOK ? inp.Parameters_OUT : _uB;
            DiffusivityA = Diffusivity(difusivityArguments_IN, _Grad_uA, inp.X);
            DiffusivityB = Diffusivity(difusivityArguments_OUT, _Grad_uB, inp.X);

            foreach (var Diffusivity in new double[]{DiffusivityA, DiffusivityB})
            {
                Debug.Assert(!double.IsNaN(Diffusivity));
                Debug.Assert(!double.IsInfinity(Diffusivity));
            }

            for (int d = 0; d < inp.D; d++) {
                // consistency term
                Acc += 0.5 * (DiffusivityA * _Grad_uA[i, d] + DiffusivityB * _Grad_uB[i, d]) * (_vA - _vB) * inp.Normal[d];
                // symmetry term                
                Acc += 0.5 * (DiffusivityA * _Grad_vA[d] + DiffusivityB * _Grad_vB[d]) * (_uA[i] - _uB[i]) * inp.Normal[d];
            }
            // penalty term          
            DiffusivityMax = (Math.Abs(DiffusivityA) > Math.Abs(DiffusivityB)) ? DiffusivityA : DiffusivityB;

            Acc -= (_uA[i] - _uB[i]) * (_vA - _vB) * pnlty * DiffusivityMax;

            return -Acc;
        }

        /// <summary>
        ///   The BoundaryEdgeForm for Dirichlet b.c. with value u_D
        /// </summary>
        protected double BoundaryEdgeFormDirichlet(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA, double u_D) {
            double Acc = 0.0;

            double pnlty = 2 * GetPenalty(inp.jCellIn, -1);
            double[] difusivityArguments_IN = prmsOK ? inp.Parameters_IN : _uA;

            double DiffusivityA = Diffusivity(difusivityArguments_IN, _Grad_uA, inp.X);
            Debug.Assert(!double.IsNaN(DiffusivityA));
            Debug.Assert(!double.IsInfinity(DiffusivityA));


            // inhom. Dirichlet b.c.
            // =====================
            for (int d = 0; d < m_D; d++) {
                Acc += (DiffusivityA * _Grad_uA[i, d]) * (_vA) * inp.Normal[d];
                Acc += (DiffusivityA * _Grad_vA[d]) * (_uA[i] - u_D) * inp.Normal[d];
            }

            Acc -= DiffusivityA * (_uA[i] - u_D) * (_vA - 0) * pnlty;
            return -Acc;
        }
 
        /// <summary>
        ///   The BoundaryEdgeForm for Neumann conditions
        /// </summary>
        protected double BoundaryEdgeFormNeumann() {
            return 0.0;
        }

        /// <summary>
        ///   BoundaryEdgeForm, depending on your physical problem you need to override this with BoundaryEdgeFormDirichlet or BoundaryEdgeFormNeumann.
        /// </summary>
        public abstract double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA);

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double Acc = 0;
            double[] difusivityArguments = prmsOK ? cpv.Parameters : U;
            double DiffusivityValue = Diffusivity(difusivityArguments, GradU, cpv.Xglobal);
            Debug.Assert(!double.IsNaN(DiffusivityValue));
            Debug.Assert(!double.IsInfinity(DiffusivityValue));

            for(int d = 0; d < cpv.D; d++)
                Acc -= DiffusivityValue * GradU[i, d] * GradV[d];
            return -Acc;
        }

        /// <summary>
        /// Arguments
        /// </summary>
        public abstract IList<string> ArgumentOrdering { get; }

        /// <summary>
        /// Parameters at linearization point to calculate material properties.
        /// </summary>
        public abstract IList<string> ParameterOrdering { get; }

    }
}
