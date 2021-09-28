using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    class Continuity : BulkEquation {

        string spcName;

        public Continuity(string spcName, int D) {
            this.spcName = spcName;
            for(int i = 0; i < D; ++i) {
                //*
                string variableName = ZwoLevelSetSolver.VariableNames.DisplacementVector(D)[i];
                AddVariableNames(variableName);
                var divergence = new Divergence(spcName, variableName, i);
                AddComponent(divergence);
                //*/

                /*
                string variableName1 = BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[i];
                AddVariableNames(variableName1);
                var divergence1 = new Divergence(spcName, variableName1, i);
                AddComponent(divergence1);
                //*/
            }
            //*
            string pressure = BoSSS.Solution.NSECommon.VariableNames.Pressure;
            AddVariableNames(pressure);
            var pressurePenalty = new EdgePenaltyForm(spcName, pressure, -1);
            AddComponent(pressurePenalty);
            //*/
        }

        public override string SpeciesName => spcName;

        public override double MassScale => 0.0;

        public override string CodomainName => BoSSS.Solution.NSECommon.EquationNames.ContinuityEquation;
    }

    class Divergence : IEdgeForm, IVolumeForm, ISpeciesFilter, ISupportsJacobianComponent {

        string speciesName;
        string[] variables;
        int d;

        public Divergence(string speciesName, string variableName, int d) {
            this.speciesName = speciesName;
            this.variables = new string[] { variableName };
            this.d = d;
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
            return flux * _vA;

        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double flux = (_uIN[0] - _uOUT[0]) * inp.Normal[d];
            flux *= 0.5 * (_vIN + _vOUT);

            //double flux = 0.5 * (_uIN[0] + _uOUT[0]) * inp.Normal[d];
            //flux *= (_vIN - _vOUT);
            return flux;

        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double flux = GradU[0,d] * V;
            return -flux;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}
