using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.MomentumEquations {
    class Gradient_pressure_xiMomentum :
            BoSSS.Foundation.IEdgeForm, // edge integrals
            BoSSS.Foundation.IVolumeForm     // volume integrals
    {
        public Gradient_pressure_xiMomentum(int _d) {
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
                //return TermActivationFlags.AllOn;
                return TermActivationFlags.UxGradV | TermActivationFlags.UxV;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                //return (TermActivationFlags.AllOn);
                return (TermActivationFlags.UxV);
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                //return TermActivationFlags.AllOn;
                return TermActivationFlags.UxV;
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

            double f_function = Globals.f_function_(r);

            double B_term;
            if (Globals.coordSys == Globals.CoordSys.hel) {
                B_term = Globals.B_term_(r);
            } else if (Globals.coordSys == Globals.CoordSys.cart) {
                B_term = 1;
            } else {
                throw new NotImplementedException("choose coordinate system!!");
            }


            double Acc = 0;

            Acc -= (-1.0 / B_term) * Psi[0] * f_function * GradV[1];  // Term 28 

            // Note: xi-derivative of f_function is 0.0!


            return (-1.0) * Acc; // Multiplikation mit -1, da die NS-Gleichung=RHS
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

            double B_term;
            if (Globals.coordSys == Globals.CoordSys.hel) {
                B_term = Globals.B_term_(r);
            } else if (Globals.coordSys == Globals.CoordSys.cart) {
                B_term = 1;
            } else {
                throw new NotImplementedException("choose coordinate system!!");
            }

            double f_function = Globals.f_function_(r);

            double Acc = 0;

            Acc += (-1.0 / B_term) * 0.5 * (Psi_IN[0] + Psi_OT[0]) * inp.Normal[d] * f_function * (V_IN - V_OT);  // Term 27

            return (-1.0) * Acc; // Multiplikation mit -1, da die NS-Gleichung=RHS
        }


        public double BoundaryEdgeForm(ref CommonParamsBnd inp,
            double[] Psi_IN, double[,] GradPsi_IN, double V_IN, double[] GradV_OT) {

            double r = inp.X[0];
            double xi = inp.X[1];

            double B_term;
            if (Globals.coordSys == Globals.CoordSys.hel) {
                B_term = Globals.B_term_(r);
            } else if (Globals.coordSys == Globals.CoordSys.cart) {
                B_term = 1;
            } else {
                throw new NotImplementedException("choose coordinate system!!");
            }
            double f_function = Globals.f_function_(r);
            double Acc = 0;

            if (Globals.AtZeroRadius(r)) {
                //
                // inner edge : multiplier f(r)=B^2
                //
                Acc += -B_term * Psi_IN[0] * inp.Normal[1] * V_IN; // Term 27
                                                                   // Nur die Innenwerte, weil Boundary

                // Wieso verwendet man nur im innneren f(r) = B(r)^2?
            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Dirichlet) {

                Acc += (-1.0 / B_term) * Psi_IN[0] * inp.Normal[1] * f_function * V_IN;   // Term 27
                                                                                          // Nur die Innenwerte, weil Boundary
            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Neumann) {
                Acc += 0;
            } else {
                throw new NotImplementedException();
            }

            return (-1.0) * Acc; // Multiplikation mit -1, da die NS-Gleichung=RHS
        }
    }
}







