﻿using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.MomentumEquations {
    class convectiveRmom : IEdgeForm, // edge integrals
                         IVolumeForm, // volume integrals
                        ISupportsJacobianComponent 
{
        public convectiveRmom() {

        }
        // Velocity components

        public IList<String> ArgumentOrdering {
            get { return new string[] { "Velocity_R", "Velocity_XI", "Velocity_ETA" }; }
        }

        public TermActivationFlags VolTerms {
            get {
                //return TermActivationFlags.AllOn;
                return TermActivationFlags.UxGradV | TermActivationFlags.UxV | TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                //return TermActivationFlags.AllOn;
                return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                //return TermActivationFlags.AllOn;
                return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | TermActivationFlags.V | TermActivationFlags.GradV;
            }
        }


        /// The parameter list:
        public IList<string> ParameterOrdering => null; 
        //{
        //    get { return new string[] { "Velocity_R0", "Velocity_XI0", "Velocity_ETA0" }; }
        //}


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double Acc = 0;
            double r = cpv.Xglobal[0];
            double a = Globals.a;
            double b = Globals.b;

            double B_term = Globals.B_term_(r);
            double f_function = Globals.f_function_(r);

            Acc += U[1] * GradU[0, 0] * f_function * V; // Term 4
            Acc += 1.0 / B_term * U[1] * GradU[0, 1] * f_function * V; // Term 5
            Acc += -1.0 * B_term * B_term / r * (b / r * U[1] + a * U[2]) * (b / r * U[1] + a * U[2]) * f_function * V;  // Term 1

            return Acc;
        }

        //public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
        //    for(int i = 0; i < 3; i++) {
        //        Console.WriteLine($"Remember: u_r Convektive is null !!!!!!");
        //    }
        //    return 0;
        //}

        public double InnerEdgeForm(ref CommonParams inp, double[] Uin, double[] Uout, double[,] _Grad_uIN, double[,] _Grad_uOUT, double V_IN, double V_OT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            // Terme 2 & 3 Wurden hier abgedeckt. Herleitung siehe unten mit Vergleich zum Paper

            double Acc = 0;
            double Flux = 0;
            double Influx = 0;
            double Outflux = 0;
            double urVel_IN = Uin[0];
            double uxiVel_IN = Uin[1];
            double urVel_OT = Uout[0];
            double uxiVel_OT = Uout[1];
            double r = inp.X[0];

            double B_term = Globals.B_term_(r);
            double f_function = Globals.f_function_(r);

            if (urVel_IN * inp.Normal[0] + uxiVel_IN * inp.Normal[1] > 0) {  // Upwind FLux!
                Flux += (urVel_IN * urVel_IN) * inp.Normal[0] + (1.0 / B_term) * (uxiVel_IN * urVel_IN) * inp.Normal[1];
                // Erster Term aus Gleichung 3.6 mit Upwind flux
                // Da Audruck positiv, deswegen Inner_Values
            } else {
                Flux += (urVel_OT * urVel_OT) * inp.Normal[0] + (1.0 / B_term) * (uxiVel_OT * urVel_OT) * inp.Normal[1];
                // Erster Term aus Gleichung 3.6 mit Upwind flux. 
                // Da Audruck negativ, deswegen Outer_Values
            }

            Influx += (urVel_IN * urVel_IN) * inp.Normal[0] + (1.0 / B_term) * (uxiVel_IN * urVel_IN) * inp.Normal[1];
            // Zweiter Term aus GLeichung 3.8
            Outflux += (urVel_OT * urVel_OT) * inp.Normal[0] + (1.0 / B_term) * (uxiVel_OT * urVel_OT) * inp.Normal[1];
            // Erster  Term aus GLeichung 3.9 bzw. Dritter Term aus Gleichung 3.8
            // Zusammenfassung der Teilterme
            // Vorzeichen nach Acc nicht ganz klar! += oder -= 
            Acc += ((Flux - Influx) * V_IN - (Flux - Outflux) * V_OT) * f_function;  // 

            
            return Acc;
        }


        double[] UDiri(double[] X) {
            double r = X[0];
            double xi = X[1];


            return new double[] {
                Globals.DirichletValue_uR(X),
                Globals.DirichletValue_uXi(X),
                Globals.DirichletValue_uEta(X)
            };


        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] GradUin, double Vin, double[] GradVin) {
            // Terme 2 & 3 Wurden hier abgedeckt. Herleitung siehe unten mit Vergleich zum Paper

            double Acc = 0;

            double[] UD;
            UD = UDiri(inp.X);
            


            double Flux = 0;
            double Influx = 0;
            double Outflux = 0;
            double urVel_IN = Uin[0];
            double uxiVel_IN = Uin[1];

            double urVel_OT = UD[0];        
            double uxiVel_OT = UD[1];       

            double r = inp.X[0];

            double B_term = Globals.B_term_(r);
            double f_function = Globals.f_function_(r);


            if (Globals.AtZeroRadius(r)) {
                //
                // inner edge
                //
                // Bei der inneren Edge für r=0 wird die Funktion f(r) durch aktiv durch B(r)^2 ersetzt,
                // um Singularitäten zu eleminieren. An den anderen Kanten wird dies nicht explizit geamcht. 
                // Wieso eigentlich nicht????

                // Term 1 Gleichung 3.6. Im Zweiten Abschnitt wurde f_function mit B_term gekürzt. f(r)=B(r)^2
                // Upwind FLux
                if (urVel_IN * inp.Normal[0] + uxiVel_IN * inp.Normal[1] > 0) { // Upwind_Flux
                    // Erster Term aus Gleichung 3.6 mit Upwind flux
                    // Da Audruck positiv, deswegen Inner_Values
                    Flux += (urVel_IN * urVel_IN) * inp.Normal[0] * f_function + B_term * (uxiVel_IN * urVel_IN) * inp.Normal[1];
                } else {
                    // Erster Term aus Gleichung 3.6 mit Upwind flux
                    // Da Audruck negativ, deswegen Outer_Values
                    Flux += (urVel_OT * urVel_OT) * inp.Normal[0] * f_function + B_term * (uxiVel_OT * urVel_OT) * inp.Normal[1];
                }

                Influx += (urVel_IN * urVel_IN) * inp.Normal[0] * f_function + B_term * (uxiVel_IN * urVel_IN) * inp.Normal[1];
                // Zweiter Term aus GLeichung 3.8
                Outflux += (urVel_OT * urVel_OT) * inp.Normal[0] * f_function + B_term * (uxiVel_OT * urVel_OT) * inp.Normal[1];
                // Dritter Term aus GLeichung 3.8 bzw. erster Term aus Fleichung 3.9
                Acc += (Flux - Influx) * Vin;
                // Boudary, deswegen nur die Innenwerte!



            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Dirichlet) {

                if (urVel_IN * inp.Normal[0] + uxiVel_IN * inp.Normal[1] > 0) {
                    // Erster Term aus Gleichung 3.6 mit Upwind flux
                    // Da Audruck positiv, deswegen Inner_Values
                    Flux += (urVel_IN * urVel_IN) * inp.Normal[0] + (1.0 / B_term) * (uxiVel_IN * urVel_IN) * inp.Normal[1];
                } else {
                    // Erster Term aus Gleichung 3.6 mit Upwind flux
                    // Da Audruck negativ, deswegen Outer_Values
                    Flux += (urVel_OT * urVel_OT) * inp.Normal[0] + (1.0 / B_term) * (uxiVel_OT * urVel_OT) * inp.Normal[1];
                }




                Influx += (urVel_IN * urVel_IN) * inp.Normal[0] + (1.0 / B_term) * (uxiVel_IN * urVel_IN) * inp.Normal[1];
                // Zweiter Term aus GLeichung 3.8
                Outflux += (urVel_OT * urVel_OT) * inp.Normal[0] + (1.0 / B_term) * (uxiVel_OT * urVel_OT) * inp.Normal[1];
                // Dritter Term aus GLeichung 3.8 bzw. erster Term aus Fleichung 3.9
                Acc += (Flux - Influx) * Vin * f_function;
                // Boudary, deswegen nur die Innenwerte!
            }

            return Acc;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var DerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var DerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { DerivEdg, DerivVol };
        }

    }

}