using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    class Continuity : BulkEquation {

        internal static bool ContinuityStabilization = false;

        string spcName;

        public Continuity(string spcName, int D, Solid Material, bool VelocityContinuity) {
            this.spcName = spcName;
            if(VelocityContinuity) {
                for(int i = 0; i < D; ++i) {
                    string velocity = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[i];
                    AddVariableNames(velocity);
                    var divergence1 = new Divergence(spcName, velocity, i);
                    AddComponent(divergence1);
                }
                //var divergence2 = new ConvectionDivergence(spcName, D);
                //AddComponent(divergence2);
            } else {
                string[] displacement = ZwoLevelSetSolver.VariableNames.DisplacementVector(D);
                for(int i = 0; i < D; ++i) {
                    AddVariableNames(displacement);
                    var divergence = new Divergence(spcName, displacement[i], i, 1.0);
                    AddComponent(divergence);
                }
                //var divergence1 = new DisplacementDivergence(spcName, displacement, 1.0);
                //AddComponent(divergence1);
            }

            if(ContinuityStabilization) {
                string pressure = BoSSS.Solution.NSECommon.VariableNames.Pressure;
                AddVariableNames(pressure);
                var pressurePenalty = new EdgePenaltyForm(spcName, pressure, - 1/ (Material.Lame2 + Material.Viscosity)); // Must scale with viscosity, see Die Pietro
                AddComponent(pressurePenalty);
            }
        }

        public override string SpeciesName => spcName;

        public override double MassScale => 0.0;

        public override string CodomainName => BoSSS.Solution.NSECommon.EquationNames.ContinuityEquation;
    }

    class Divergence : IEdgeForm, IVolumeForm, ISpeciesFilter, ISupportsJacobianComponent {

        string speciesName;
        string[] variables;
        int d;
        double scale;
        
        public Divergence(string speciesName, string variableName, int d, double scale = 1.0) {
            this.speciesName = speciesName;
            this.variables = new string[] { variableName };
            this.d = d;
            this.scale = scale;
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public IList<string> ArgumentOrdering => variables;

        public IList<string> ParameterOrdering => null;

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.GradUxV; }
        }

        public string ValidSpecies => speciesName;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double flux = _uA[0] * inp.Normal[d];
            return scale * flux * _vA;
            //return 0.0;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double flux = (_uIN[0] - _uOUT[0]) * inp.Normal[d];
            flux *=  0.5 * (_vIN + _vOUT);

            //double flux = 0.5 * (_uIN[0] + _uOUT[0]) * inp.Normal[d];
            //flux *= (_vIN - _vOUT);
            return scale * flux;

        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double flux = GradU[0,d] * V;
            return - scale * flux;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    class ConvectionDivergence : IEdgeForm, IVolumeForm, ISpeciesFilter, ISupportsJacobianComponent {

        string speciesName;
        string[] variables;
        int D;

        public ConvectionDivergence(string speciesName, int D) {
            this.speciesName = speciesName;
            this.variables = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D).Cat(VariableNames.DisplacementVector(D));
            this.D = D;
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.None; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.GradUxV; }
        }

        public IList<string> ArgumentOrdering => variables;

        public IList<string> ParameterOrdering => null;

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.UxGradV|TermActivationFlags.GradUxGradV; }
        }

        public string ValidSpecies => speciesName;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0.0;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc = 0;
            for(int i = 0; i < D; i++) {
                double t = 0;
                for(int j = 0; j < D; ++j) {
                    t = _uIN[j]  * (_Grad_uIN[D+i, j]) + (_uOUT[j] * _Grad_uOUT[D+i, j]);
                }
                acc += 0.5 * t * inp.Normal[i];
            }
            return -acc * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int i = 0; i < D; i++) {
                double t = 0;
                for(int j = 0; j < D; ++j) {
                    t = U[j] * GradU[D + i, j];
                }
                acc += t * GradV[i];
            }
            return acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DivergenceDerivEdg, DivergenceDerivVol };
        }
    }
}
