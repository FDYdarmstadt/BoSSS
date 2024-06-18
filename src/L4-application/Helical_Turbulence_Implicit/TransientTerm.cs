using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit {
    class TransientTerm : IVolumeForm {


        double factor;
        int velIndex;

        public TransientTerm(double _factor, int _velIndex) {
            this.factor = _factor;
            this.velIndex = _velIndex;
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { (new string[] { "Velocity_R", "Velocity_XI", "Velocity_ETA" })[velIndex] }; }
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[0];
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double r = cpv.Xglobal[0];
            return U[0] * V * factor * Globals.f_function_(r);

        }
    }
}
