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
    class SurfacePenaltyForm : IEdgeForm {

        string[] variableNames;
        double penaltyScale;

        public SurfacePenaltyForm(string variableName, double penaltyScale) {
            this.variableNames = new string[] { variableName };
            this.penaltyScale = penaltyScale;
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => null;

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.None;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, 
            double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            double acc = penaltyScale * (_uIN[0] - _uOUT[0]) * (_vIN - _vOUT);
            return acc;
        }
    }
}