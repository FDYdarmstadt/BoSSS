using BoSSS.Platform.LinAlg;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE {
    public static class Globals {
        public static double Det2x2(this double[,] M) {
            return M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0];
        }

        public static int GetBDFOrder(this HelicalControl c) {
            int bdfOrder;

            if(c.TimeSteppingScheme == TimeSteppingScheme.CrankNicolson) {
                bdfOrder = -1;
            } else if(c.TimeSteppingScheme == TimeSteppingScheme.ImplicitEuler) {
                bdfOrder = 1;
            } else if(c.TimeSteppingScheme.ToString().StartsWith("BDF")) {
                bdfOrder = Convert.ToInt32(c.TimeSteppingScheme.ToString().Substring(3));
            } else {
                throw new NotImplementedException("todo");
            }
            if(c.BurstSaves < bdfOrder) {
                c.BurstSaves= bdfOrder;
            }
            return bdfOrder;
        }

        public static double a = 1.0;
        public const double b = 1.0;
        public const double nu = 1.0;

        public const double G = 1.0;
        public const double K = 1.0;
        public const double Z = 1.0;

        /// <summary>
        /// parameters for Backward-Difference Method for time stepping
        /// </summary>
        public const double beta0 = 11.0 / 6.0;
        public const double beta1 = 3.0;
        public const double beta2 = -3.0 / 2.0;
        public const double beta3 = 1.0 / 3.0;

        public const double g1 = 3;
        public const double g2 = -3;
        public const double g3 = 1;

        public static double gammaUR { get; internal set; } // hier bessere aber komplizierte variante penalty term
        public static double gammaUXI { get; internal set; }
        public static double gammaETA { get; internal set; }

        static public Func<double[], BoundaryTypeE> BoundaryType;
        static public Func<double[], double> DirichletValue_uR;
        static public Func<double[], double> DirichletValue_uXi;
        static public Func<double[], double> DirichletValue_uEta;
        static public Func<double[], double> etaMomSourceTerm;


        public static double B_term_(double r) {
            //double last_r = -1;
            //double last_form = double.NaN;
            double last_form;

            //if (r == last_r) {
            //    return last_form;
            //} else {
            //    last_r = r;
            if(b == 0 && r == 0) {
                last_form = 1.0 / a;
            } else {
                last_form = r / Math.Sqrt(a * a * r * r + b * b);
            }

            return last_form;
            //}
        }

        public static bool AtZeroRadius(double r) {
            if (r < -1.0e-9)
                throw new ArgumentException("Grid at negative radius - cant be!");
            return (r.Abs() < 1.0e-9);
        }

        public static CoordSys coordSys = CoordSys.hel;

        public enum CoordSys {
            pol,
            hel,
            cart
        }
        public enum Multiplier {
            one,            // takes the constant f(r)=1 as the multipier
            Bsq             // takes f(r)=B^2 as the multiplier
        }

        public static Multiplier activeMult {
            get { 
                return m_activeMult; 
            }
            set {
                m_activeMult = value;
            }
        }

        static Multiplier m_activeMult = Multiplier.one;

        public static double f_function_(double r) {
            double f_function;
            switch (activeMult) {
                case Multiplier.one:
                    f_function = 1.0;
                    return f_function;
                case Multiplier.Bsq:
                    f_function = Globals.B_term_(r) * Globals.B_term_(r);
                    return f_function;
                default:
                    throw new NotImplementedException("missing multiplier " + activeMult);
            }

        }

        public static double df_function_(double r) {
            double df_function;
            switch (activeMult) {
                case Multiplier.one:
                    df_function = 0.0;
                    return df_function;
                case Multiplier.Bsq:
                    df_function = 2 * r / (a * a * r * r + b * b) - 2 * r * r * r * a * a / ((a * a * r * r + b * b) * (a * a * r * r + b * b));
                    return df_function;
                default:
                    throw new NotImplementedException("missing derivative of multiplier ");
            }

        }
        public static double penaltyScaling(double r) {
            switch(activeMult) {
                case Multiplier.one:
                return 1;
                case Multiplier.Bsq:
                return r;
                default:
                throw new NotImplementedException("missing derivative of multiplier ");
            }
        }
        public static bool steady = false;
        public static bool pressureStabilConti = true;
        public static bool pressureStabilEtaMom = false;

        public static double MaxAmp = 5.0;
    }


}
