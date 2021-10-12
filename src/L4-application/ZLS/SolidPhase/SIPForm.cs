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
    /// For debugging only, to be delted
    /// </summary>
    class Fake_ipFlux : BoSSS.Solution.NSECommon.SIPLaplace, ISpeciesFilter {

       public Fake_ipFlux(double penalty_const, string varName, double visc)
            : base(penalty_const, varName) //
        {
            m_visc = visc;
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

        //public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.None;

        //public override TermActivationFlags InnerEdgeTerms => base.InnerEdgeTerms;

        //public override TermActivationFlags VolTerms => base.VolTerms;

        
    }


    /// <summary>
    /// Viscosity/energy dissipation in the solid phase
    /// </summary>
    public class SIPForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent, IEquationComponentCoefficient {
        double viscosity;
        string species;
        int d;
        string[] variableNames;
        public double PenaltySafety = 1.3;

        public SIPForm(string species, string[] variables, int d, double viscosity) {
            this.species = species;
            this.viscosity = viscosity;
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
