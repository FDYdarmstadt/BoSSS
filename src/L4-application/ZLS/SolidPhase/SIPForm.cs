using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {

    /// <summary>
    /// Viscosity/energy dissipation in the solid phase
    /// </summary>
    public class SIPForm : IVolumeForm, IEdgeForm, ISpeciesFilter, ISupportsJacobianComponent, IEquationComponentCoefficient {
        protected double viscosity;
        string species;
        protected int d;
        protected string[] variableNames;
        IncompressibleMultiphaseBoundaryCondMap boundaryMap;
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

        public SIPForm(string species, string[] variables, int d, double viscosity, IncompressibleMultiphaseBoundaryCondMap boundaryMap, double __PenaltySafety = 4.0) 
            :this(species, variables, d, viscosity, __PenaltySafety) {
            this.boundaryMap = boundaryMap;
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxGradV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => new string[] { };

        public TermActivationFlags BoundaryEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.V | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | TermActivationFlags.V | TermActivationFlags.GradV); }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return (TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV); }
        }

        public string ValidSpecies => species;

        public virtual double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN) {
            if(boundaryMap != null) {
                IncompressibleBcType edgType = boundaryMap.EdgeTag2Type[inp.EdgeTag];
                double acc1 = 0.0;
                double pnlty = PenaltyIn(inp.jCellIn);
                Vector dirichlet = new Vector(D);
                for(int i = 0; i < D; ++i) {
                    dirichlet[i] = boundaryMap.bndFunction[variableNames[i]][inp.EdgeTag](inp.X, inp.time);
                }
                switch(edgType) {
                    case IncompressibleBcType.Wall:
                    for(int i = 0; i < D; i++) {
                        acc1 -= viscosity * _Grad_uIN[d, i] * _vIN * inp.Normal[i];  // consistency term  
                        acc1 -= viscosity * _Grad_vIN[i] * (_uIN[d] - dirichlet[d]) * inp.Normal[i];  // symmetry term
                    }
                    acc1 += PenaltySafety * pnlty* (_uIN[d] - dirichlet[d]) * _vIN * viscosity;
                    break;
                    case IncompressibleBcType.FreeSlip:
                    for(int i = 0; i < D; i++) {
                        for(int j = 0; j < D; ++j) {
                            acc1 -= viscosity * inp.Normal[i] * _Grad_uIN[i, j] * inp.Normal[j] * _vIN * inp.Normal[d];  // consistency term  
                            acc1 -= viscosity * inp.Normal[d] * _Grad_vIN[j] * inp.Normal[j] * (_uIN[i] - dirichlet[i]) * inp.Normal[i];  // symmetry term
                        }
                        acc1 += PenaltySafety * pnlty * (_uIN[i] - dirichlet[i]) * inp.Normal[i] * _vIN * inp.Normal[d];
                    }
                    break;
                    case IncompressibleBcType.Outflow:
                    case IncompressibleBcType.Pressure_Outlet:
                    case IncompressibleBcType.Velocity_Inlet:
                    break;
                    default:
                    throw new NotImplementedException();
                }
                return acc1;
            } else {
                Vector dirichlet = new Vector(D);
                double acc1 = 0.0;
                double pnlty = PenaltyIn(inp.jCellIn);
                for(int i = 0; i < D; i++) {
                    acc1 -= viscosity * _Grad_uIN[d, i] * _vIN * inp.Normal[i];  // consistency term  
                    acc1 -= viscosity * _Grad_vIN[i] * (_uIN[d] - dirichlet[d]) * inp.Normal[i];  // symmetry term
                }
                acc1 += PenaltySafety * (_uIN[d] - dirichlet[d]) * _vIN * pnlty * viscosity;
                return acc1;
            }
        }

        MultidimensionalArray cj;

        double penalty;

        protected int D;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();
            D = cs.GrdDat.SpatialDimension;
            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes
            penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            this.cj = cs.CellLengthScales;
            //this.cj = ((BoSSS.Foundation.Grid.Classic.GridData)(cs.GrdDat)).Cells.CellLengthScale;
        }

        protected double PenaltyIn(int jCellIn) {
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
