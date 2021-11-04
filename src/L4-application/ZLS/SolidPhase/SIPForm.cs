using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {

    /// <summary>
    /// Symmetric interior penalty with only penalties
    /// </summary>
    class Penalty_ipFlux : BoSSS.Solution.NSECommon.SIPLaplace, ISpeciesFilter {

       public Penalty_ipFlux(double penalty_const, string varName, double visc)
            : base(penalty_const, varName) //
        {
            m_visc = visc;
            base.m_alpha = 0.0;
        }

     
        public string ValidSpecies => "C";

        protected override double g_Diri(ref CommonParamsBnd inp) {
            return 0;
        }

        protected override double g_Neum(ref CommonParamsBnd inp) {
            return 0;
        }

        double m_visc;

        public override double Nu(double[] x, double[] p, int jCell) {
            return -m_visc;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            return true;
        }

        /*
        override public double InnerEdgeForm(ref BoSSS.Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
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


            //for(int d = 0; d < inp.D; d++) {
            //    Acc -= 1.0 * (nuA * _Grad_uA[0, d] - nuB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normal[d];  // consistency term
            //    //Acc += 0.5 * (nuA * _Grad_vA[d] + nuB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normal[d];  // symmetry term
            //}



            return Acc;

        }

        /// <summary>
        /// Integrand on boundary mesh edges of the SIP
        /// </summary>
        override public double BoundaryEdgeForm(ref BoSSS.Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
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

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.AllOn;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.AllOn;
        */
    }





    /// <summary>
    /// Viscosity/energy dissipation in the solid phase
    /// </summary>
    public class SIPForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent, IEquationComponentCoefficient {
        double viscosity;
        string species;
        int d;
        string[] variableNames;
        public double PenaltySafety;

        public SIPForm(string species, string[] variables, int d, double viscosity, double __PenaltySafety = 4.0) {
            this.species = species;
            this.viscosity = viscosity;
            this.PenaltySafety = __PenaltySafety;
            if(this.viscosity <= 0.0) {
                throw new ArgumentException($"Viscosity must be positive, but got a value of {this.viscosity}");
            }
            this.variableNames = variables;
            this.d = d;
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxGradV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[] { };

        public TermActivationFlags BoundaryEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | TermActivationFlags.V | TermActivationFlags.GradV); }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV); }
        }

        public string ValidSpecies => species;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            double acc1 = 0.0;
            double pnlty = PenaltyIn(inp.jCellIn);

            Vector dirichlet = new Vector(D);

            for (int i = 0; i < D; i++) {
                acc1 -= viscosity * _Grad_uIN[d, i] * _vIN  * inp.Normal[i];  // consistency term  
                acc1 -= viscosity * _Grad_vIN[i] * (_uIN[d] - dirichlet[d]) * inp.Normal[i];  // symmetry term
            }
            acc1 += PenaltySafety * (_uIN[d] - dirichlet[d]) * _vIN  * pnlty * viscosity;
            return acc1;
        }

        MultidimensionalArray cj;

        double penalty;

        int D;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
            D = cs.GrdDat.SpatialDimension;
            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes
            penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            this.cj = cs.CellLengthScales;
            //this.cj = ((BoSSS.Foundation.Grid.Classic.GridData)(cs.GrdDat)).Cells.CellLengthScale;
        }

        double PenaltyIn(int jCellIn) {
            double penaltySizeFactor_A = 1.0 / cj[jCellIn];
            
            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor_A * penalty;
            if(µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }

        double PenaltyOut(int jCellOut) {
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0;
            
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penalty));

            double µ = penaltySizeFactor_B * penalty;
            if (µ.IsNaNorInf())
                throw new ArithmeticException("Inf/NaN in penalty computation.");
            return µ;
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc1 = 0.0;
            for(int i = 0; i < D; i++) {
                acc1 -= 0.5 * viscosity * ( _Grad_uIN[d, i] + _Grad_uOUT[d, i]) * (_vIN - _vOUT) * inp.Normal[i];  // consistency term  
                acc1 -= 0.5 * viscosity * (_Grad_vIN[i] + _Grad_vOUT[i]) * (_uIN[d] - _uOUT[d]) * inp.Normal[i];  // symmetry term
            }
            double penalty = Math.Max(PenaltyIn(inp.jCellIn), PenaltyOut(inp.jCellOut));
            acc1 += (_uIN[d] - _uOUT[d]) * PenaltySafety * penalty *(_vIN - _vOUT) * viscosity;
            return acc1;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int i = 0; i < D; ++i) {
                acc += GradU[d, i] * GradV[i];
            }
            acc *= viscosity;
            return acc ;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}
