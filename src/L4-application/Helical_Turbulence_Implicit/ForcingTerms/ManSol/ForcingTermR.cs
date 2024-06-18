using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.ForcingTerms.ManSol {
    class ForcingTermR : IVolumeForm {

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
                    return -1 * (0.2e1 * r * ((nu * (a * a * r * r + b * b) * Math.Sin(xi) + (-0.2e1 * Math.Pow(r, 0.3e1) + 0.2e1 * r) * Math.Pow(Math.Cos(xi), 0.2e1)) * Math.Exp(-r * r) + (0.2e1 * Math.Pow(r, 0.3e1) - r) * Math.Pow(Math.Cos(xi), 0.2e1) * Math.Exp(-0.2e1 * r * r) - nu * (a * a * r * r + b * b) * Math.Sin(xi) - r * Math.Pow(Math.Cos(xi), 0.2e1)) * a * b * Math.Sqrt(a * a * r * r + b * b) + (-(-0.4e1 * a * a * nu * Math.Pow(r, 0.6e1) - 0.2e1 * a * a * Math.Pow(r, 0.5e1) + nu * (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(r, 0.4e1) - 0.2e1 * b * b * Math.Pow(r, 0.3e1) + 0.2e1 * ((a * a + 0.4e1) * b * b + a * a / 0.2e1) * nu * r * r + Math.Pow(b, 0.4e1) * nu - b * b * nu) * (a * a * r * r + b * b) * Math.Sin(xi) + 0.2e1 * r * (-b * b * ((a * a + 0.2e1) * r * r + b * b - 0.1e1) * Math.Pow(Math.Cos(xi), 0.2e1) + Math.Pow(a * a * r * r + b * b, 0.2e1) * r * r)) * Math.Exp(-r * r) - 0.2e1 * r * (-(-0.4e1 * Math.Pow(r, 0.4e1) + (a * a + 0.4e1) * r * r + b * b - 0.1e1) * b * b * Math.Pow(Math.Cos(xi), 0.2e1) / 0.2e1 + Math.Pow(a * a * r * r + b * b, 0.2e1) * r * r) * Math.Exp(-0.2e1 * r * r) + (Math.Pow(a, 0.4e1) * Math.Pow(r, 0.4e1) + (0.2e1 * a * a * b * b + a * a) * r * r + Math.Pow(b, 0.4e1) - b * b) * (a * a * r * r + b * b) * nu * Math.Sin(xi) + Math.Pow(Math.Cos(xi), 0.2e1) * b * b * r * (a * a * r * r + b * b - 0.1e1)) * Math.Pow(a * a * r * r + b * b, -0.2e1) * Math.Pow(r, -0.2e1) * V;
                    case Globals.Multiplier.Bsq:
                    return -1 * (0.2e1 * a * ((nu * (a * a * r * r + b * b) * Math.Sin(xi) + (-0.2e1 * Math.Pow(r, 0.3e1) + 0.2e1 * r) * Math.Pow(Math.Cos(xi), 0.2e1)) * Math.Exp(-r * r) + (0.2e1 * Math.Pow(r, 0.3e1) - r) * Math.Pow(Math.Cos(xi), 0.2e1) * Math.Exp(-0.2e1 * r * r) - nu * (a * a * r * r + b * b) * Math.Sin(xi) - r * Math.Pow(Math.Cos(xi), 0.2e1)) * b * r * Math.Sqrt(a * a * r * r + b * b) + (-(a * a * r * r + b * b) * (-0.4e1 * a * a * nu * Math.Pow(r, 0.6e1) - 0.2e1 * a * a * Math.Pow(r, 0.5e1) + nu * (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(r, 0.4e1) - 0.2e1 * b * b * Math.Pow(r, 0.3e1) + 0.2e1 * ((a * a + 0.4e1) * b * b + a * a / 0.2e1) * nu * r * r + Math.Pow(b, 0.4e1) * nu - b * b * nu) * Math.Sin(xi) + 0.2e1 * (-((a * a + 0.2e1) * r * r + b * b - 0.1e1) * b * b * Math.Pow(Math.Cos(xi), 0.2e1) + Math.Pow(a * a * r * r + b * b, 0.2e1) * r * r) * r) * Math.Exp(-r * r) - 0.2e1 * (-b * b * (-0.4e1 * Math.Pow(r, 0.4e1) + (a * a + 0.4e1) * r * r + b * b - 0.1e1) * Math.Pow(Math.Cos(xi), 0.2e1) / 0.2e1 + Math.Pow(a * a * r * r + b * b, 0.2e1) * r * r) * r * Math.Exp(-0.2e1 * r * r) + nu * (a * a * r * r + b * b) * (Math.Pow(a, 0.4e1) * Math.Pow(r, 0.4e1) + (0.2e1 * a * a * b * b + a * a) * r * r + Math.Pow(b, 0.4e1) - b * b) * Math.Sin(xi) + Math.Pow(Math.Cos(xi), 0.2e1) * b * b * r * (a * a * r * r + b * b - 0.1e1)) * Math.Pow(a * a * r * r + b * b, -0.3e1) * V;
                    default:
                    throw new NotImplementedException("missing multiplier" + Globals.activeMult);
                }
            } else {
                switch(Globals.activeMult) {
                    case Globals.Multiplier.one:
                    // Copy DDD
                    return -1 * (0.2e1 * r * Math.Cos(t) * (((-0.2e1 * Math.Pow(r, 0.3e1) + 0.2e1 * r) * Math.Pow(Math.Cos(xi), 0.2e1) * Math.Cos(t) + Math.Sin(xi) * nu * (a * a * r * r + b * b)) * Math.Exp(-r * r) + 0.2e1 * r * (r * r - 0.1e1 / 0.2e1) * Math.Cos(t) * Math.Pow(Math.Cos(xi), 0.2e1) * Math.Exp(-0.2e1 * r * r) - r * Math.Cos(t) * Math.Pow(Math.Cos(xi), 0.2e1) - Math.Sin(xi) * nu * (a * a * r * r + b * b)) * b * a * Math.Sqrt(a * a * r * r + b * b) + (0.2e1 * (-((a * a + 0.2e1) * r * r + b * b - 0.1e1) * b * b * Math.Pow(Math.Cos(xi), 0.2e1) + Math.Pow(a * a * r * r + b * b, 0.2e1) * r * r) * r * Math.Pow(Math.Cos(t), 0.2e1) - (a * a * r * r + b * b) * Math.Sin(xi) * (-0.4e1 * a * a * nu * Math.Pow(r, 0.6e1) - 0.2e1 * a * a * Math.Pow(r, 0.5e1) + nu * (Math.Pow(a, 0.4e1) + 0.4e1 * a * a - 0.4e1 * b * b) * Math.Pow(r, 0.4e1) - 0.2e1 * b * b * Math.Pow(r, 0.3e1) + 0.2e1 * nu * ((a * a + 0.4e1) * b * b + a * a / 0.2e1) * r * r + Math.Pow(b, 0.4e1) * nu - b * b * nu) * Math.Cos(t) + Math.Sin(xi) * r * r * Math.Sin(t) * Math.Pow(a * a * r * r + b * b, 0.2e1)) * Math.Exp(-r * r) - 0.2e1 * r * Math.Pow(Math.Cos(t), 0.2e1) * (-(-0.4e1 * Math.Pow(r, 0.4e1) + (a * a + 0.4e1) * r * r + b * b - 0.1e1) * b * b * Math.Pow(Math.Cos(xi), 0.2e1) / 0.2e1 + Math.Pow(a * a * r * r + b * b, 0.2e1) * r * r) * Math.Exp(-0.2e1 * r * r) + b * b * r * Math.Pow(Math.Cos(xi), 0.2e1) * (a * a * r * r + b * b - 0.1e1) * Math.Pow(Math.Cos(t), 0.2e1) + (a * a * r * r + b * b) * Math.Sin(xi) * nu * (Math.Pow(a, 0.4e1) * Math.Pow(r, 0.4e1) + (0.2e1 * a * a * b * b + a * a) * r * r + Math.Pow(b, 0.4e1) - b * b) * Math.Cos(t) - Math.Sin(xi) * r * r * Math.Sin(t) * Math.Pow(a * a * r * r + b * b, 0.2e1)) * Math.Pow(a * a * r * r + b * b, -0.2e1) * Math.Pow(r, -0.2e1) * V;
                    case Globals.Multiplier.Bsq:
                    return -1 * (2 * b * nu * r * (a * Math.Cos(t) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) - b * Math.Exp(-r * r) * Math.Cos(t) * Math.Sin(xi) * (Math.Exp(r * r) + 2 * r * r - 1) / (r * Math.Pow(a * a * r * r + b * b, 1 / 2)))) * Math.Pow(a * a * r * r + b * b, -(3 / 2)) - nu * Math.Cos(t) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) - r * r * r * Math.Pow(a * Math.Cos(t) * Math.Cos(xi) * (Math.Exp(-r * r) - 1) - b * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(xi) * (Math.Exp(r * r) + 2 * r * r - 1) * Math.Pow(r * Math.Pow(a * a * r * r + b * b, 1 / 2), -1), 2) * Math.Pow(a * a * r * r + b * b, -2) + r * r * Math.Sin(t) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) * Math.Pow(a * a * r * r + b * b, -1) - nu * Math.Cos(t) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) * Math.Pow(a * a * r * r + b * b, -1) + 2 * r * r * r * Math.Exp(-r * r) * Math.Cos(t) * Math.Sin(xi) * Math.Pow(a * a * r * r + b * b, -1) - 2 * r * r * r * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(t) * Math.Sin(xi) * Math.Sin(xi) * (Math.Exp(-r * r) - 1) / (a * a * r * r + b * b) - r * Math.Exp(-r * r) * Math.Cos(t) * Math.Cos(t) * Math.Cos(xi) * Math.Cos(xi) * (Math.Exp(-r * r) - 1) * (Math.Exp(r * r) + 2 * r * r - 1) * Math.Pow(a * a * r * r + b * b, -1) + 4 * nu * r * r * Math.Exp(-r * r) * Math.Cos(t) * Math.Sin(xi) * (r * r - 1) * Math.Pow(a * a * r * r + b * b, -1) * V;
                    default:
                    throw new NotImplementedException("missing multiplier" + Globals.activeMult);
                }
            }

            //return -ExactResidual(cpv.Xglobal, t) * V;


        }
    }
}
