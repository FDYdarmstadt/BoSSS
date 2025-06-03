using BoSSS.Foundation;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace StokesHelical_Ak.ForcingTerms.Hagen_Poiseulle {
    class ForcingTermXi_Hagen : IVolumeForm {

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


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double t = cpv.time;
            double a = Globals.a;
            double b = Globals.b;
            double r = cpv.Xglobal[0];
            double xi = cpv.Xglobal[1];
            double nu = Globals.nu;
            double maxAmp = Globals.MaxAmp;

            return maxAmp * (r * a * a / Math.Sqrt(a * a * r * r + b * b)) * Globals.f_function_(r) * V;
        }
    }
}
