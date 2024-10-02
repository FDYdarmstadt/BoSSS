using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.ContinuityEquation {
    // ###################################### helical formulation ########################################### //

    // U[0] = ur
    // U[1] = uxi
    // U[2] = ueta


    public class Conti :
            BoSSS.Foundation.IEdgeForm, // edge integrals
            BoSSS.Foundation.IVolumeForm,     // volume integrals
            BoSSS.Foundation.ISupportsJacobianComponent // For Jacobian, required for Newton Solver

    {

        public Conti(int _ResolutionXI) {
            this.ResolutionXI = _ResolutionXI;
        }

        /// The resolution in r-direction: needed for pressure stabilization
        int ResolutionXI;

        /// The parameter list for the divergence is empty:
        public IList<string> ParameterOrdering {
            get { return null; }
        }

        public IList<String> ArgumentOrdering {
            get { return new string[] { "Velocity_R", "Velocity_XI", "Pressure" }; }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxGradV | TermActivationFlags.UxV;
                //return TermActivationFlags.AllOn;
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                //return TermActivationFlags.AllOn;
                return TermActivationFlags.UxV;
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
                //return TermActivationFlags.AllOn;
            }
        }



        public double VolumeForm(ref CommonParamsVol cpv,
            double[] U, double[,] GradU,
            double V, double[] GradV) {

            double Acc = 0;

            double r = cpv.Xglobal[0];
            double xi = cpv.Xglobal[1];
            double B_term;
            if (Globals.coordSys == Globals.CoordSys.hel) {
                B_term = Globals.B_term_(r);
            } else if (Globals.coordSys == Globals.CoordSys.cart) {
                B_term = 1;
            } else {
                throw new NotImplementedException("choose coordinate system!!");
            };


            double f_function = Globals.f_function_(r);
            double df_function = Globals.df_function_(r);



            Acc -= U[0] * (df_function * V + f_function * GradV[0]);    // Term 3 und 2   // term diff(UR(r,xi),r)
            if (Globals.coordSys == Globals.CoordSys.hel) {
                Acc += 1 / r * U[0] * f_function * V;                   // Term 1         // source term 1/r*UR(r,xi)
            }
            Acc -= U[1] * f_function / B_term * GradV[1];               // Term 4         // term diff(UXI(r,xi),xi)



            return Acc;
        }



        public double InnerEdgeForm(ref CommonParams inp,
            double[] U_IN, double[] U_OT, double[,] GradU_IN, double[,] GradU_OT,
            double V_IN, double V_OT, double[] GradV_IN, double[] GradV_OT) {

            double Acc = 0;

            double r = inp.X[0];

            double B_term;
            if (Globals.coordSys == Globals.CoordSys.hel) {
                B_term = Globals.B_term_(r);
            } else if (Globals.coordSys == Globals.CoordSys.cart) {
                B_term = 1;
            } else {
                throw new NotImplementedException("choose coordinate system!!");
            }
            double f_function = Globals.f_function_(r);


            Acc += 0.5 * (U_IN[0] + U_OT[0]) * f_function * inp.Normal[0] * (V_IN - V_OT);            // Term 5      // term diff(UR(r,xi),r) 
            Acc += 0.5 * (U_IN[1] + U_OT[1]) * f_function / B_term * inp.Normal[1] * (V_IN - V_OT);   // Term 6      // term diff(UXI(r,xi),xi)


            // pressure stabilization term, U[2] is the pressure
            //Acc += 1 / Globals.nu * 2 * Math.PI / 16 * (U_IN[2] - U_OT[2]) * (V_IN - V_OT) * Globals.penaltyScaling(r);
            //Acc += 1 / Globals.nu * 2 * Math.PI / ResolutionXI * (U_IN[2] - U_OT[2]) * (V_IN - V_OT) * Globals.penaltyScaling(r);
              

            if(Globals.pressureStabilConti == true) {
                double h_edge = (inp.GridDat as GridData).Edges.h_max_Edge[inp.iEdge];

                // Have to be true. See Details on DiPietro 2013#
                double factor = Globals.penaltyScaling(r);
                if(Globals.activeMult == Globals.Multiplier.Bsq) {
                    factor *= r;
                }
                Acc += 1 / Globals.nu * 2 * Math.PI / ResolutionXI * (U_IN[2] - U_OT[2]) * (V_IN - V_OT) * factor;
            }



            return Acc;
        }



        double[] UDiri(double[] X) {
            double r = X[0];
            double xi = X[1];

            return new double[] {
                Globals.DirichletValue_uR(X),
                Globals.DirichletValue_uXi(X),
            };

        }



        public double BoundaryEdgeForm(ref CommonParamsBnd inp,
        double[] U_IN, double[,] GradU_IN, double V_IN, double[] GradV_IN) {

            double Acc = 0;
            double r = inp.X[0];
            double[] UD;
            UD = UDiri(inp.X);


            double B_term;
            if (Globals.coordSys == Globals.CoordSys.hel) {
                B_term = Globals.B_term_(r);
            } else if (Globals.coordSys == Globals.CoordSys.cart) {
                B_term = 1;
            } else {
                throw new NotImplementedException("choose coordinate system!!");
            }
            double f_function = Globals.f_function_(r);

            if (Globals.AtZeroRadius(r)) {
                //
                // inner edge
                //
                Acc += U_IN[0] * f_function * inp.Normal[0] * V_IN;         // Term 5  // term diff(UR(r,xi),r) 


                //Acc += U_IN[1] * f_function / B_term * inp.Normale[1] * V_IN;        // term diff(UXI(r,xi),xi)
                Acc += U_IN[1] * B_term * inp.Normal[1] * V_IN;             // Term 6  // term diff(UXI(r,xi),xi)         // multiplicator f(r)=B^2


                // All Dead Code! Never true. See Details on DePietro 2013
                //if(Globals.pressureStabilConti == true) {
                //    Acc += 1 / Globals.nu * 2 * Math.PI / ResolutionXI * (U_IN[2] - UD[2]) * (V_IN) * Globals.penaltyScaling(r);
                //}

            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Dirichlet) {

                Acc += U_IN[0] * f_function * inp.Normal[0] * V_IN;           // Term 5     // term diff(UR(r,xi),r) 
                Acc += U_IN[1] * f_function / B_term * inp.Normal[1] * V_IN;  // Term 6     // term diff(UXI(r,xi),xi)


                //if(Globals.pressureStabilConti == true) {
                //    // All Dead Code! Never true. See Details on DePietro 2013
                //    Acc += 1 / Globals.nu * 2 * Math.PI / ResolutionXI * (U_IN[2] - UD[2]) * (V_IN) * Globals.penaltyScaling(r);
                //}


            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Neumann) {
                Acc += 0;
            } else {
                throw new NotImplementedException("Wrong boundary type");
            }

            return Acc;
        }


        /// <summary>
        /// Linear component - derivative is just this.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new[] { this };
        }

    }
}
