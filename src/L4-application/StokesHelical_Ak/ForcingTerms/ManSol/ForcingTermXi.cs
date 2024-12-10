using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace StokesHelical_Ak.ForcingTerms.ManSol {
    class ForcingTermXi : IVolumeForm {

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
                switch(Globals.activeMult) {
                    case Globals.Multiplier.one:
                    return -1 * 0.2e1 * Math.Pow(a * a * r * r + b * b, -0.5e1 / 0.2e1) * Math.Cos(xi) * (-0.2e1 * (((a * a * r * r + b * b) * Math.Sin(xi) + r * nu * (a * a * r * r + a * a / 0.2e1 + b * b)) * Math.Exp(-r * r) - (a * a * r * r + b * b) * Math.Sin(xi) * Math.Exp(-0.2e1 * r * r) / 0.2e1 + Math.Sin(xi) * (-a * a * r * r / 0.2e1 - b * b / 0.2e1) - a * a * r * nu / 0.2e1) * r * r * a * b * Math.Sqrt(a * a * r * r + b * b) + (-0.2e1 * Math.Pow(r, 0.3e1) * (r - 0.1e1) * (r + 0.1e1) * (a * a * r * r + a * a + b * b) * (a * a * r * r + b * b) * Math.Sin(xi) - 0.4e1 * Math.Pow(a, 0.4e1) * nu * Math.Pow(r, 0.10e2) + a * a * nu * (Math.Pow(a, 0.4e1) + 0.10e2 * a * a - 0.8e1 * b * b) * Math.Pow(r, 0.8e1) - Math.Pow(a, 0.6e1) * Math.Pow(r, 0.7e1) / 0.2e1 - (Math.Pow(a, 0.6e1) + (-0.6e1 * b * b + 0.2e1) * Math.Pow(a, 0.4e1) - 0.48e2 * a * a * b * b + 0.8e1 * Math.Pow(b, 0.4e1)) * nu * Math.Pow(r, 0.6e1) / 0.2e1 - 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(r, 0.5e1) - ((b * b - 0.1e1) * Math.Pow(a, 0.4e1) + (-0.6e1 * Math.Pow(b, 0.4e1) + 0.4e1 * b * b) * a * a - 0.28e2 * Math.Pow(b, 0.4e1)) * nu * Math.Pow(r, 0.4e1) / 0.2e1 - 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(r, 0.3e1) + ((b * b - 0.4e1) * a * a + 0.2e1 * Math.Pow(b, 0.4e1) - 0.10e2 * b * b) * b * b * nu * r * r / 0.2e1 - Math.Pow(b, 0.6e1) * r / 0.2e1 + Math.Pow(b, 0.6e1) * nu / 0.2e1 - Math.Pow(b, 0.4e1) * nu / 0.2e1) * Math.Exp(-r * r) - Math.Pow(r, 0.3e1) * Math.Sin(xi) * (a * a + 0.2e1 * b * b) * (a * a * r * r + b * b) * Math.Exp(-0.2e1 * r * r) + (-Math.Pow(a, 0.4e1) * Math.Pow(r, 0.5e1) - a * a * b * b * Math.Pow(r, 0.3e1)) * Math.Sin(xi) + Math.Pow(a, 0.6e1) * Math.Pow(r, 0.7e1) / 0.2e1 + Math.Pow(a, 0.6e1) * nu * Math.Pow(r, 0.6e1) / 0.2e1 + 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(r, 0.5e1) + Math.Pow(a, 0.4e1) * nu * (b - 0.1e1) * (b + 0.1e1) * Math.Pow(r, 0.4e1) / 0.2e1 + 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(r, 0.3e1) + (-Math.Pow(b, 0.4e1) * nu / 0.2e1 + 0.2e1 * b * b * nu) * a * a * r * r + Math.Pow(b, 0.6e1) * r / 0.2e1 - Math.Pow(b, 0.6e1) * nu / 0.2e1 + Math.Pow(b, 0.4e1) * nu / 0.2e1) * Math.Pow(r, -0.2e1) * V;
                    case Globals.Multiplier.Bsq:
                    return -1 * 0.2e1 * Math.Pow(a * a * r * r + b * b, -0.7e1 / 0.2e1) * Math.Cos(xi) * (-0.2e1 * a * (((a * a * r * r + b * b) * Math.Sin(xi) + r * nu * (a * a * r * r + a * a / 0.2e1 + b * b)) * Math.Exp(-r * r) - (a * a * r * r + b * b) * Math.Sin(xi) * Math.Exp(-0.2e1 * r * r) / 0.2e1 + Math.Sin(xi) * (-a * a * r * r / 0.2e1 - b * b / 0.2e1) - a * a * r * nu / 0.2e1) * b * r * r * Math.Sqrt(a * a * r * r + b * b) + (-0.2e1 * Math.Pow(r, 0.3e1) * (r - 0.1e1) * (r + 0.1e1) * (a * a * r * r + a * a + b * b) * (a * a * r * r + b * b) * Math.Sin(xi) - 0.4e1 * Math.Pow(a, 0.4e1) * nu * Math.Pow(r, 0.10e2) + a * a * nu * (Math.Pow(a, 0.4e1) + 0.10e2 * a * a - 0.8e1 * b * b) * Math.Pow(r, 0.8e1) - Math.Pow(a, 0.6e1) * Math.Pow(r, 0.7e1) / 0.2e1 - nu * (Math.Pow(a, 0.6e1) + (-0.6e1 * b * b + 0.2e1) * Math.Pow(a, 0.4e1) - 0.48e2 * a * a * b * b + 0.8e1 * Math.Pow(b, 0.4e1)) * Math.Pow(r, 0.6e1) / 0.2e1 - 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(r, 0.5e1) - nu * ((b * b - 0.1e1) * Math.Pow(a, 0.4e1) + (-0.6e1 * Math.Pow(b, 0.4e1) + 0.4e1 * b * b) * a * a - 0.28e2 * Math.Pow(b, 0.4e1)) * Math.Pow(r, 0.4e1) / 0.2e1 - 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(r, 0.3e1) + ((b * b - 0.4e1) * a * a + 0.2e1 * Math.Pow(b, 0.4e1) - 0.10e2 * b * b) * nu * b * b * r * r / 0.2e1 - Math.Pow(b, 0.6e1) * r / 0.2e1 + Math.Pow(b, 0.6e1) * nu / 0.2e1 - Math.Pow(b, 0.4e1) * nu / 0.2e1) * Math.Exp(-r * r) - Math.Pow(r, 0.3e1) * Math.Sin(xi) * (a * a + 0.2e1 * b * b) * (a * a * r * r + b * b) * Math.Exp(-0.2e1 * r * r) + (-Math.Pow(a, 0.4e1) * Math.Pow(r, 0.5e1) - a * a * b * b * Math.Pow(r, 0.3e1)) * Math.Sin(xi) + Math.Pow(a, 0.6e1) * Math.Pow(r, 0.7e1) / 0.2e1 + Math.Pow(a, 0.6e1) * nu * Math.Pow(r, 0.6e1) / 0.2e1 + 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(r, 0.5e1) + Math.Pow(a, 0.4e1) * nu * (b - 0.1e1) * (b + 0.1e1) * Math.Pow(r, 0.4e1) / 0.2e1 + 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(r, 0.3e1) + (0.2e1 * b * b * nu - Math.Pow(b, 0.4e1) * nu / 0.2e1) * a * a * r * r + Math.Pow(b, 0.6e1) * r / 0.2e1 - Math.Pow(b, 0.6e1) * nu / 0.2e1 + Math.Pow(b, 0.4e1) * nu / 0.2e1) * V;
                    default:
                    throw new NotImplementedException("missing multiplier" + Globals.activeMult);
                }
            } else {
                switch(Globals.activeMult) {
                    case Globals.Multiplier.one:
                    // Copy DDD
                    return -1 * (0.2e1 * (r * r * ((-0.2e1 * Math.Sin(xi) * (a * a * r * r + b * b) * Math.Cos(t) - 0.2e1 * r * nu * (a * a * r * r + a * a / 0.2e1 + b * b)) * Math.Exp(-r * r) + Math.Sin(xi) * (a * a * r * r + b * b) * Math.Cos(t) * Math.Exp(-0.2e1 * r * r) + Math.Sin(xi) * (a * a * r * r + b * b) * Math.Cos(t) + a * a * r * nu) * Math.Cos(t) * b * a * Math.Sqrt(a * a * r * r + b * b) + (-0.2e1 * Math.Sin(xi) * Math.Pow(r, 0.3e1) * (r - 0.1e1) * (r + 0.1e1) * (a * a * r * r + b * b) * (a * a * r * r + a * a + b * b) * Math.Pow(Math.Cos(t), 0.2e1) + (-0.4e1 * Math.Pow(a, 0.4e1) * nu * Math.Pow(r, 0.10e2) + a * a * nu * (Math.Pow(a, 0.4e1) + 0.10e2 * a * a - 0.8e1 * b * b) * Math.Pow(r, 0.8e1) - Math.Pow(a, 0.6e1) * Math.Pow(r, 0.7e1) / 0.2e1 - (Math.Pow(a, 0.6e1) + (-0.6e1 * b * b + 0.2e1) * Math.Pow(a, 0.4e1) - 0.48e2 * a * a * b * b + 0.8e1 * Math.Pow(b, 0.4e1)) * nu * Math.Pow(r, 0.6e1) / 0.2e1 - 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(r, 0.5e1) - nu * ((b * b - 0.1e1) * Math.Pow(a, 0.4e1) + (-0.6e1 * Math.Pow(b, 0.4e1) + 0.4e1 * b * b) * a * a - 0.28e2 * Math.Pow(b, 0.4e1)) * Math.Pow(r, 0.4e1) / 0.2e1 - 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(r, 0.3e1) + ((b * b - 0.4e1) * a * a + 0.2e1 * Math.Pow(b, 0.4e1) - 0.10e2 * b * b) * nu * b * b * r * r / 0.2e1 - Math.Pow(b, 0.6e1) * r / 0.2e1 + Math.Pow(b, 0.6e1) * nu / 0.2e1 - Math.Pow(b, 0.4e1) * nu / 0.2e1) * Math.Cos(t) - Math.Pow(a * a * r * r + b * b, 0.2e1) * Math.Sin(t) * r * r * (r * r - 0.1e1 / 0.2e1)) * Math.Exp(-r * r) - Math.Sin(xi) * Math.Pow(r, 0.3e1) * Math.Pow(Math.Cos(t), 0.2e1) * (a * a + 0.2e1 * b * b) * (a * a * r * r + b * b) * Math.Exp(-0.2e1 * r * r) - a * a * Math.Sin(xi) * Math.Pow(r, 0.3e1) * (a * a * r * r + b * b) * Math.Pow(Math.Cos(t), 0.2e1) + (Math.Pow(a, 0.6e1) * Math.Pow(r, 0.7e1) / 0.2e1 + Math.Pow(a, 0.6e1) * nu * Math.Pow(r, 0.6e1) / 0.2e1 + 0.3e1 / 0.2e1 * Math.Pow(a, 0.4e1) * b * b * Math.Pow(r, 0.5e1) + Math.Pow(a, 0.4e1) * nu * (b - 0.1e1) * (b + 0.1e1) * Math.Pow(r, 0.4e1) / 0.2e1 + 0.3e1 / 0.2e1 * a * a * Math.Pow(b, 0.4e1) * Math.Pow(r, 0.3e1) + (-Math.Pow(b, 0.4e1) * nu / 0.2e1 + 0.2e1 * b * b * nu) * a * a * r * r + Math.Pow(b, 0.6e1) * r / 0.2e1 - Math.Pow(b, 0.6e1) * nu / 0.2e1 + Math.Pow(b, 0.4e1) * nu / 0.2e1) * Math.Cos(t) - r * r * Math.Sin(t) * Math.Pow(a * a * r * r + b * b, 0.2e1) / 0.2e1) * Math.Cos(xi) * Math.Pow(a * a * r * r + b * b, -0.5e1 / 0.2e1) * Math.Pow(r, -0.2e1)) * V;
                    case Globals.Multiplier.Bsq:
                    return -1 * (r * r * Math.Sin(t) * (Math.Cos(xi) * (Math.Exp(-r * r) - 1) * Math.Pow(a * a * r * r + b * b, -(1 / 2)) - 2 * r * r * Math.Exp(-r * r) * Math.Cos(xi) * Math.Pow(a * a * r * r + b * b, -(1 / 2)))) * Math.Pow(a * a * r * r + b * b, -1) - 2 * b * nu * r * r * (a * a * a * r * Math.Cos(t) * Math.Cos(xi) * (Math.Exp(-r * r) - 1) * Math.Pow(a * a * r * r + b * b, -(3 / 2)) - b * Math.Cos(t) * Math.Cos(xi) * (Math.Exp(-r * r) - 1) * Math.Pow(r * r, -1) + 2 * a * r * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(xi) * Math.Pow(a * a * r * r + b * b, -(1 / 2))) * Math.Pow(a * a * r * r + b * b, -(3 / 2)) - r * Math.Cos(t) * Math.Cos(xi) * (Math.Exp(-r * r) - 1) * Math.Pow(a * a * r * r + b * b, -(1 / 2)) + nu * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(xi) * (Math.Exp(r * r) + 2 * r * r - 1) * Math.Pow(a * a * r * r + b * b, -(1 / 2)) - nu * r * r * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(xi) * (12 * b * b * b * b + 2 * a * a * b * b - a * a * a * a * r * r + 2 * a * a * a * a * r * r * r * r - 20 * a * a * a * a * r * r * r * r * r * r + 8 * a * a * a * a * r * r * r * r * r * r * r * r - 28 * b * b * b * b * r * r + 8 * b * b * b * b * r * r * r * r - 2 * a * a * b * b * Math.Exp(r * r) + a * a * a * a * r * r * Math.Exp(r * r) + 8 * a * a * b * b * r * r - 48 * a * a * b * b * r * r * r * r + 16 * a * a * b * b * r * r * r * r * r * r) * Math.Pow(a * a * r * r + b * b, -(7 / 2)) - r * Math.Exp(-2 * r * r) * Math.Cos(t) * Math.Cos(t) * Math.Cos(xi) * Math.Sin(xi) * (Math.Exp(r * r) + 2 * r * r - 1) * (Math.Exp(r * r) + 2 * r * r - 1) * Math.Pow(a * a * r * r + b * b, -(3 / 2)) + b * b * nu * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(xi) * (2 * a * a * r * r + b * b) * (Math.Exp(r * r) + 2 * r * r - 1) * Math.Pow(a * a * r * r + b * b, -(7 / 2)) + 2 * a * b * r * r * Math.Cos(t) * Math.Cos(t) * Math.Cos(xi) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) * (Math.Exp(-r * r) - 1) * Math.Pow(a * a * r * r + b * b, -2) + r * r * r * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(t) * Math.Cos(xi) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) * (a * a * Math.Exp(r * r) - a * a - 6 * b * b - 4 * a * a * r * r + 4 * a * a * r * r * r * r + 4 * b * b * r * r) * Math.Pow(a * a * r * r + b * b, -(5 / 2)) - b * b * r * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(t) * Math.Cos(xi) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) * (Math.Exp(r * r) + 2 * r * r - 1) * Math.Pow(a * a * r * r + b * b, -(5 / 2)) * V;
                    default:
                    throw new NotImplementedException("missing multiplier" + Globals.activeMult);
                }
            }

            //return -ExactResidual(cpv.Xglobal, t) * V;


        }
    }
}
