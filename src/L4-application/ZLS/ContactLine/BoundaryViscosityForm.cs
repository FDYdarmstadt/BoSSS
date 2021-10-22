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
    class BoundaryViscosityForm : IVolumeForm, IEdgeForm {

        string[] variableNames;
        int D;
        int d;
        double viscosity;

        public BoundaryViscosityForm(string[] variableNames, int d, int D, double viscosity) {
            this.variableNames = variableNames;
            this.D = D;
            this.viscosity = viscosity;
            this.d = d;
        }

        public TermActivationFlags VolTerms => TermActivationFlags.AllOn;

        public IList<string> ArgumentOrdering => variableNames;

        public IList<string> ParameterOrdering {
            get {
                string[] parameters = NSEVariableNames.AsLevelSetVariable(
                    VariableNames.SolidLevelSetCG, NSEVariableNames.NormalVector(D)).ToArray();
                return parameters;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.None;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.AllOn;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            
            Vector Nsurf_IN = new Vector(inp.Parameters_IN).Normalize();
            Vector Nsurf_OUT = new Vector(inp.Parameters_OUT).Normalize();

            double[,] Psurf_IN = SurfaceUtilities.SurfaceProjection(Nsurf_IN);
            double[,] Psurf_OUT = SurfaceUtilities.SurfaceProjection(Nsurf_OUT);

            Vector tauL_IN = SurfaceUtilities.Tangent(Nsurf_IN, inp.Normal);
            Vector tauL_OUT = SurfaceUtilities.Tangent(Nsurf_OUT, inp.Normal);

            double divUsurf_IN = 0.0;
            double divUsurf_OUT = 0.0;
            for (int i = 0; i < inp.D; i++) {
                for (int j = 0; j < inp.D; j++) {
                    divUsurf_IN += Psurf_IN[i, j] * _Grad_uIN[j, i];
                    divUsurf_OUT += Psurf_OUT[i, j] * _Grad_uOUT[j, i];
                }
            }

            double acc = 0.0;
            acc -= 0.5 * viscosity * (divUsurf_IN * tauL_IN[d] + divUsurf_OUT * tauL_OUT[d]) * (_vIN - _vOUT);
            
            acc += 100 * viscosity * (_uIN[d] - _uOUT[d]) * (_vIN - _vOUT);
            return acc;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            int D = cpv.D;

            Vector flux = new Vector(D);
            
            Vector Nsurf = new Vector(cpv.Parameters).Normalize();

            double[,] Psurf = SurfaceUtilities.SurfaceProjection(Nsurf);

            double surfDiv = 0.0;
            for (int i = 0; i < D; i++) {
                for (int j = 0; j < D; j++) {
                    surfDiv += Psurf[i, j] * GradU[j, i];
                }
            }

            for (int i = 0; i < D; i++) {
                flux[i] = viscosity * surfDiv * Psurf[d, i];
            }

            double acc = 0;
            for (int i = 0; i < D; i++)
                acc += flux[i] * GradV[i];
            
            //acc = 0;
            return acc;
        }
    }
}
