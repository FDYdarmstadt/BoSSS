using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.MomentumEquations {
    class Gradient_pressure_rMomentum :
            BoSSS.Foundation.IEdgeForm, // edge integrals
            BoSSS.Foundation.IVolumeForm     // volume integrals
    {
        public Gradient_pressure_rMomentum(int _d) {
            this.d = _d;
        }

        /// The component index of the gradient:
        int d;

        /// As ususal, we do not use parameters:
        public IList<string> ParameterOrdering {
            get { return null; }
        }

        /// We have one argument, the pressure $\psi$:
        public IList<String> ArgumentOrdering {
            get { return new string[] { "Pressure" }; }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.AllOn;
                //return TermActivationFlags.UxGradV | TermActivationFlags.UxV;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.AllOn);
                //return (TermActivationFlags.UxV);
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.AllOn;
                //return TermActivationFlags.UxV;
            }
        }

        /// The volume integrand, for a vector-valued test-function $\vec{v}$
        /// would be $-\operatorname{div}{\vec{v}} \psi$. Our test function $v$
        /// is scalar-valued, so e.g. for $\code{d} = 0$ we have
        /// $\vec{v} = (v,0)$. In this case, our volume integrand reduces as 
        /// $-\operatorname{div}{\vec{v}} \psi = -\partial_x v \psi$:
        public double VolumeForm(ref CommonParamsVol cpv,
               double[] Psi, double[,] GradPsi,
               double V, double[] GradV) {

            double r = cpv.Xglobal[0];
            double xi = cpv.Xglobal[1];

            double df_function = Globals.df_function_(r);
            double f_function = Globals.f_function_(r);

            double Acc = 0;

            Acc -= -Psi[0] * (df_function * V + f_function * GradV[d]); // Term 25 und 26

            return (-1.0) * Acc;  // Mulitplikation mit -1, da die NS-Gleichung=RHS!
        }

        /// On interior cell edges, we simply use a central-difference flux.
        /// Again, we consider a scalar test function, so we have
        /// $ \jump{\psi} \vec{v} \cdot \vec{n} = \jump{\psi} v n_d $,
        /// where $n_d$ is the $d$--th component of $\vec{n}$:
        public double InnerEdgeForm(ref CommonParams inp,
            double[] Psi_IN, double[] Psi_OT,
            double[,] GradPsi_IN, double[,] GradPsi_OT,
            double V_IN, double V_OT, double[] GradV_IN, double[] GradV_OT) {
            double r = inp.X[0];
            double xi = inp.X[1];

            double f_function = Globals.f_function_(r);
            double Acc = 0;

            Acc += -0.5 * (Psi_IN[0] + Psi_OT[0]) * inp.Normal[0] * f_function * (V_IN - V_OT); // Term 24



            return (-1.0) * Acc; // Multiplikation mit -1, da die NS-Gleichung=RHS
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp,
            double[] Psi_IN, double[,] GradPsi_IN, double V_IN, double[] GradV_OT) {

            double r = inp.X[0];
            double xi = inp.X[1];
            double f_function = Globals.f_function_(r);



            double Acc = 0;

            if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Dirichlet) {

                Acc += -Psi_IN[0] * inp.Normal[0] * f_function * V_IN  ;   //  Term 24
                // Boundary, deswegen nur die Innenwerte
                // Wieso gibt es hier keine Implementierung für die Werte and er Inneren Kante? 
                // Nicht nötig, es gibt kein 1/B(r). Deswegen kann man f(r) = B(r)^2 nicht verwenden....
            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Neumann) {
                Acc += 0;
            } else {
                throw new NotImplementedException();
            }

            return (-1.0) * Acc;  // Multiplikation mit -1, da die NS-Gleichung=RHS

        }

    }
}