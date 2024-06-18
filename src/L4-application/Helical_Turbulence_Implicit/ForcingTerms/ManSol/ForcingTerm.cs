using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.ForcingTerms.ManSol {
    class ForcingTerm : IVolumeForm {

        public IList<string> ArgumentOrdering {
            get {
                return new string[0];
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[0];
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }


        public Func<double[], double, double> ExactResidual;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double t = cpv.time;
            return -ExactResidual(cpv.Xglobal, t) * V;
        }
    }
}
