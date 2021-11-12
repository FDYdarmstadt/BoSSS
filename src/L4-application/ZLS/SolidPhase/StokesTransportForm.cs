using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {

    //All boundaries are a Wall for now
    public class StokesTransportForm : IVolumeForm, IEdgeForm, ISpeciesFilter {

        string speciesName;

        double rho;

        string[] variableNames;

        string[] parameternames;

        int D;

        int d;

        public StokesTransportForm(string speciesName, string[] variableNames, int d, int D, double rho) {
            this.speciesName = speciesName;
            this.variableNames = variableNames;
            this.D = D;
            this.d = d;
            this.rho = rho;
            this.parameternames = ZwoLevelSetSolver.VariableNames.Displacement0Vector(D);
        }

        public TermActivationFlags VolTerms {
            get { return TermActivationFlags.UxGradV | TermActivationFlags.GradV; }
        }

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering => parameternames;

        public TermActivationFlags BoundaryEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public TermActivationFlags InnerEdgeTerms {
            get { return TermActivationFlags.UxV | TermActivationFlags.V; }
        }

        public string ValidSpecies => speciesName;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uA, double _vIN, double[] _Grad_vA) {

            Vector VelocityIn = new Vector(inp.Parameters_IN, 0, D);

            // Upwinding:
            if (VelocityIn * inp.Normal >= 0)
            {
                return rho * _uIN[d] * (VelocityIn * inp.Normal) * (_vIN); // outflow
            }
            else
            {
                return 0.0 * (_vIN); // inflow
            }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double r = 0.0;

            Vector VelocityIn = new Vector(inp.Parameters_IN, 0, D);
            Vector VelocityOt = new Vector(inp.Parameters_OUT, 0, D);
            Vector VelocityAvg = 0.5 * (VelocityIn + VelocityOt);

            // Upwinding:
            if (VelocityAvg * inp.Normal >= 0)
            {
                return rho * _uIN[d] * (VelocityIn * inp.Normal) * (_vIN - _vOUT);
            }
            else
            {
                return rho * _uOUT[d] * (VelocityOt * inp.Normal) * (_vIN - _vOUT);
            }
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int i = 0; i < D; i++)
                acc += cpv.Parameters[i] * GradV[i];
            acc *= U[d];
            acc *= rho;
            return -acc;
        }
    }
}
