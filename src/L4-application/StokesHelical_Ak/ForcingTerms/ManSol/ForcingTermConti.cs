using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace StokesHelical_Ak.ForcingTerms.ManSol {
    class ForcingTermConti : IVolumeForm {

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


        //public Func<double[], double, double> ExactResidual;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double t = cpv.time;
            double a = Globals.a;
            double b = Globals.b;
            double r = cpv.Xglobal[0];
            double xi = cpv.Xglobal[1];
            double nu = Globals.nu;
            return -1 * 0 * V;

            //return -ExactResidual(cpv.Xglobal, t) * V;


        }
    }
}
