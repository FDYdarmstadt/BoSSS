using BoSSS.Foundation;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NSEVariableNames = BoSSS.Solution.NSECommon.VariableNames;

namespace ZwoLevelSetSolver.ContactLine {
    class SurfaceGradientPenaltyForm : IEdgeForm {

        string[] variableNames;
        double scale;
        int D;

        public SurfaceGradientPenaltyForm(string variableName, int D, double scale) {
            this.variableNames = new string[] { variableName };
            this.scale = scale;
            this.D = D;
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering {
            get {
                string[] parameters = NSEVariableNames.AsLevelSetVariable(
                    VariableNames.SolidLevelSetCG, NSEVariableNames.NormalVector(D)).ToArray();
                return parameters;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.None;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.GradUxV;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, 
            double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            
            Vector Nsurf_IN = new Vector(inp.Parameters_IN).Normalize();
            Vector Nsurf_OUT = new Vector(inp.Parameters_OUT).Normalize();

            Vector tauL_IN = SurfaceUtilities.Rotate(Nsurf_IN);
            Vector tauL_OUT = SurfaceUtilities.Rotate(Nsurf_OUT);

            Vector tau = 0.5 * (tauL_IN + tauL_OUT);

            double gradUsurf_IN = 0.0;
            double gradUsurf_OUT = 0.0;
            for (int i = 0; i < inp.D; i++) {
                gradUsurf_IN += tau[i] * _Grad_uIN[0,i];
                gradUsurf_OUT += tau[i] * _Grad_uOUT[0, i];
            }

            
            double acc = scale * (gradUsurf_IN - gradUsurf_OUT) * (_vIN - _vOUT);
            return acc;
        }
    }
}
