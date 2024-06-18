using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.ForcingTerms.ManSol {
    class ForcingTermEta : IVolumeForm {

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
            if(Globals.steady == true) {
                switch(Globals.activeMult) { // Navier Stokes case
                    case Globals.Multiplier.one:
                    return -1 * (-Math.Pow(a * a * r * r + b * b, -0.5e1 / 0.2e1) * Math.Cos(xi) * ((((-0.2e1 * a * a * b * b * Math.Pow(r, 0.3e1) - 0.2e1 * Math.Pow(b, 0.4e1) * r) * Math.Sin(xi) + (-0.4e1 * Math.Pow(a, 0.4e1) * Math.Pow(r, 0.8e1) + (Math.Pow(a, 0.6e1) + 0.4e1 * Math.Pow(a, 0.4e1) - 0.8e1 * a * a * b * b) * Math.Pow(r, 0.6e1) + (-0.4e1 * Math.Pow(b, 0.4e1) + (0.3e1 * Math.Pow(a, 0.4e1) + 0.8e1 * a * a) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(r, 0.4e1) + ((0.3e1 * a * a + 0.4e1) * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * r * r + Math.Pow(b, 0.6e1)) * nu) * Math.Exp(-r * r) + b * b * r * Math.Sin(xi) * (a * a * r * r + b * b) * Math.Exp(-0.2e1 * r * r) + (a * a * b * b * Math.Pow(r, 0.3e1) + Math.Pow(b, 0.4e1) * r) * Math.Sin(xi) - (Math.Pow(a, 0.6e1) * Math.Pow(r, 0.6e1) + (0.3e1 * Math.Pow(a, 0.4e1) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(r, 0.4e1) + (0.3e1 * a * a * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * r * r + Math.Pow(b, 0.6e1)) * nu) * Math.Sqrt(a * a * r * r + b * b) - 0.2e1 * r * ((-0.4e1 * a * a * Math.Pow(r, 0.6e1) + (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(r, 0.4e1) + ((0.2e1 * a * a + 0.8e1) * b * b + a * a) * r * r + Math.Pow(b, 0.4e1) - b * b) * Math.Exp(-r * r) - Math.Pow(a, 0.4e1) * Math.Pow(r, 0.4e1) + (-0.2e1 * a * a * b * b - a * a) * r * r - Math.Pow(b, 0.4e1) + b * b) * a * b * nu) * Math.Pow(r, -0.2e1)) * V;
                    case Globals.Multiplier.Bsq:
                    return -1 * (-Math.Pow(a * a * r * r + b * b, -0.7e1 / 0.2e1) * Math.Cos(xi) * ((((-0.2e1 * a * a * b * b * Math.Pow(r, 0.3e1) - 0.2e1 * Math.Pow(b, 0.4e1) * r) * Math.Sin(xi) + (-0.4e1 * Math.Pow(a, 0.4e1) * Math.Pow(r, 0.8e1) + (Math.Pow(a, 0.6e1) + 0.4e1 * Math.Pow(a, 0.4e1) - 0.8e1 * a * a * b * b) * Math.Pow(r, 0.6e1) + (-0.4e1 * Math.Pow(b, 0.4e1) + (0.3e1 * Math.Pow(a, 0.4e1) + 0.8e1 * a * a) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(r, 0.4e1) + ((0.3e1 * a * a + 0.4e1) * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * r * r + Math.Pow(b, 0.6e1)) * nu) * Math.Exp(-r * r) + b * b * r * Math.Sin(xi) * (a * a * r * r + b * b) * Math.Exp(-0.2e1 * r * r) + (a * a * b * b * Math.Pow(r, 0.3e1) + Math.Pow(b, 0.4e1) * r) * Math.Sin(xi) - nu * (Math.Pow(a, 0.6e1) * Math.Pow(r, 0.6e1) + (0.3e1 * Math.Pow(a, 0.4e1) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(r, 0.4e1) + (0.3e1 * a * a * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * r * r + Math.Pow(b, 0.6e1))) * Math.Sqrt(a * a * r * r + b * b) - 0.2e1 * a * ((-0.4e1 * a * a * Math.Pow(r, 0.6e1) + (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(r, 0.4e1) + ((0.2e1 * a * a + 0.8e1) * b * b + a * a) * r * r + Math.Pow(b, 0.4e1) - b * b) * Math.Exp(-r * r) - Math.Pow(a, 0.4e1) * Math.Pow(r, 0.4e1) + (-0.2e1 * a * a * b * b - a * a) * r * r - Math.Pow(b, 0.4e1) + b * b) * nu * b * r)) * V;
                    default:
                    throw new NotImplementedException("missing multiplier" + Globals.activeMult);
                }
            } else {
                switch(Globals.activeMult) {
                    case Globals.Multiplier.one:
                    // Copy DDD
                    return -1 * (-(((-0.2e1 * Math.Sin(xi) * b * b * r * (a * a * r * r + b * b) * Math.Pow(Math.Cos(t), 0.2e1) + (-0.4e1 * Math.Pow(a, 0.4e1) * Math.Pow(r, 0.8e1) + (Math.Pow(a, 0.6e1) + 0.4e1 * Math.Pow(a, 0.4e1) - 0.8e1 * a * a * b * b) * Math.Pow(r, 0.6e1) + (-0.4e1 * Math.Pow(b, 0.4e1) + (0.3e1 * Math.Pow(a, 0.4e1) + 0.8e1 * a * a) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(r, 0.4e1) + ((0.3e1 * a * a + 0.4e1) * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * r * r + Math.Pow(b, 0.6e1)) * nu * Math.Cos(t) - r * r * Math.Sin(t) * Math.Pow(a * a * r * r + b * b, 0.2e1)) * Math.Exp(-r * r) + Math.Sin(xi) * b * b * r * (a * a * r * r + b * b) * Math.Pow(Math.Cos(t), 0.2e1) * Math.Exp(-0.2e1 * r * r) + Math.Sin(xi) * b * b * r * (a * a * r * r + b * b) * Math.Pow(Math.Cos(t), 0.2e1) - (Math.Pow(a, 0.6e1) * Math.Pow(r, 0.6e1) + (0.3e1 * Math.Pow(a, 0.4e1) * b * b + Math.Pow(a, 0.4e1)) * Math.Pow(r, 0.4e1) + (0.3e1 * a * a * Math.Pow(b, 0.4e1) + 0.2e1 * a * a * b * b) * r * r + Math.Pow(b, 0.6e1)) * nu * Math.Cos(t) + r * r * Math.Sin(t) * Math.Pow(a * a * r * r + b * b, 0.2e1)) * Math.Sqrt(a * a * r * r + b * b) - 0.2e1 * r * nu * ((-0.4e1 * a * a * Math.Pow(r, 0.6e1) + (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(r, 0.4e1) + ((0.2e1 * a * a + 0.8e1) * b * b + a * a) * r * r + Math.Pow(b, 0.4e1) - b * b) * Math.Exp(-r * r) - Math.Pow(a, 0.4e1) * Math.Pow(r, 0.4e1) + (-0.2e1 * a * a * b * b - a * a) * r * r - Math.Pow(b, 0.4e1) + b * b) * Math.Cos(t) * b * a) * Math.Cos(xi) * Math.Pow(a * a * r * r + b * b, -0.5e1 / 0.2e1) * Math.Pow(r, -0.2e1)) * V;
                    case Globals.Multiplier.Bsq:
                    return -1 * (r * r * Math.Cos(xi) * Math.Sin(t) * (Math.Exp(-r * r) - 1)) * Math.Pow(a * a * r * r + b * b, -1) - nu * Math.Cos(t) * Math.Cos(xi) * (Math.Exp(-r * r) - 1) - 2 * r * r * r * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(t) * Math.Cos(xi) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) * Math.Pow(a * a * r * r + b * b, -1) + 4 * nu * r * r * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(xi) * (r * r - 1) * Math.Pow(a * a * r * r + b * b, -1) + a * a * r * r * r * Math.Cos(t) * Math.Cos(t) * Math.Cos(xi) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) * (Math.Exp(-r * r) - 1) * Math.Pow(a * a * r * r + b * b, -2) + r * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(t) * Math.Cos(xi) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) * (Math.Exp(r * r) + 2 * r * r - 1) * Math.Pow(a * a * r * r + b * b, -1) - a * a * nu * r * r * Math.Cos(t) * Math.Cos(xi) * (a * a * r * r + 2 * b * b) * (Math.Exp(-r * r) - 1) * Math.Pow(a * a * r * r + b * b, -3) + 2 * a * b * nu * r * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(xi) * (b * b * Math.Exp(r * r) - b * b * b * b * Math.Exp(r * r) - b * b + b * b * b * b + a * a * r * r + 4 * a * a * r * r * r * r - 4 * a * a * r * r * r * r * r * r + a * a * a * a * r * r * r * r + 8 * b * b * r * r - 4 * b * b * r * r * r * r - a * a * r * r * Math.Exp(r * r) - a * a * a * a * r * r * r * r * Math.Exp(r * r) + 2 * a * a * b * b * r * r - 2 * a * a * b * b * r * r * Math.Exp(r * r)) * Math.Pow(a * a * r * r + b * b, -(7 / 2)) * V;
                    default:
                    throw new NotImplementedException("missing multiplier" + Globals.activeMult);
                }
            }

            //return -ExactResidual(cpv.Xglobal, t) * V;


        }
    }
}
