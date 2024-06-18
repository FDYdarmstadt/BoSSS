using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IncompressibleNSE.Helical_Turbulence_Implicit.MomentumEquations {
    class xiMomentum :
            IEdgeForm, // edge integrals
            IVolumeForm     // volume integrals

    {



        public TermSwitch m_TermSwitch;



        private Func<double, int, int, MultidimensionalArray, double> m_ComputePenalty;

        private double m_penaltyFactor;


        public xiMomentum(TermSwitch termSwitch, double PenaltyFactor, Func<double, int, int, MultidimensionalArray, double> ComputePenalty = null) {
            m_TermSwitch = termSwitch;
            m_ComputePenalty = ComputePenalty;
            m_penaltyFactor = PenaltyFactor;
        }


        /// <summary>
        /// coefficient of durdxi
        /// </summary>
        /// <param name="r"></param>
        /// <returns></returns>
        static double CoeffUR(double r) {
            return 2.0 * Globals.b * Globals.b * Globals.B_term_(r) / (r * r * r);
        }

        /// <summary>
        /// coefficient of durdxi at the boundary r=0: this is coefficient of durdxi * f_function(r)
        /// </summary>
        /// <param name="r"></param>
        /// <returns></returns>
        static double CoeffURbnd(double r) {
            double a = Globals.a;
            double b = Globals.b;
            return 2.0 * Globals.b * Globals.b / (Math.Sqrt(a * a * r * r + b * b) * (a * a * r * r + b * b));
        }


        /// <summary>
        /// coefficient of duetadr
        /// </summary>
        /// <param name="r"></param>
        /// <returns></returns>
        static double CoeffUETA(double r) {
            return 2.0 * Globals.b * Globals.B_term_(r) / r;
        }

        /// <summary>
        /// coefficient of duetadr at the boundary r=0: this is coefficient of duretadr * f_function(r)
        /// </summary>
        /// <param name="r"></param>
        /// <returns></returns>
        static double CoeffUETAbnd(double r) {
            double f_function = Globals.f_function_(r);
            double a = Globals.a;
            double b = Globals.b;
            return 2.0 * Globals.b / (Math.Sqrt(a * a * r * r + b * b)) * f_function;
        }



        /// <summary>
        /// First derivative of coefficient of duetadr
        /// </summary>
        /// <param name="r"></param>
        /// <returns></returns>
        static double CoeffUETA_deriv(double r) {
            double B_term = Globals.B_term_(r);
            double a = Globals.a;
            double b = Globals.b;
            return -2.0 * b * a * a * B_term * B_term * B_term / (r * r);
        }


        /// The parameter list for the divergence is empty:
        public IList<string> ParameterOrdering {
            get { return null; }
        }


        // Velocity components

        public IList<String> ArgumentOrdering {
            get { return new string[] {"Velocity_R", "Velocity_XI", "Velocity_ETA" }; }
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


        // Volume Form


        double[,] CoeffGradUXI(double r) {

            double f_function = Globals.f_function_(r);
            double nu = Globals.nu;

            switch (Globals.coordSys) {
                case Globals.CoordSys.pol:
                    //return new double[,] { { r, 0 }, { 0, 1 / r } };
                    return new double[,] { { f_function, 0 }, { 0, f_function / (r * r) } };
                case Globals.CoordSys.hel:
                    return new double[,] { { nu * f_function, 0 }, { 0, nu * f_function / (Globals.B_term_(r) * Globals.B_term_(r)) } };
                case Globals.CoordSys.cart:
                    return new double[,] { { 1, 0 }, { 0, 1 } };
                default:
                    throw new NotImplementedException("missing coordinate system " + Globals.coordSys);
            }

        }



        public double VolumeForm(ref CommonParamsVol cpv,
        double[] U, double[,] GradU,
        double V, double[] GradV) {

            double Acc = 0;

            double r = cpv.Xglobal[0];
            double xi = cpv.Xglobal[1];


            double f_function = Globals.f_function_(r);
            double df_function = Globals.df_function_(r);

            double B_term = Globals.B_term_(r);
            double a = Globals.a;
            double nu = Globals.nu;



            if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomXI) != 0) {
                Acc -= nu * U[1] * ((df_function * r - f_function) / (r * r) * V + f_function / r * GradV[0]);
                //Acc -= nu * U[1] * df_function * r * V / (r * r); Term 18
                //Acc -= nu * U[1] * f_function * V / (r * r);      Term 19
                //Acc -= nu * U[1] * f_function * GradV[0] / r;     Term 20

                Acc += nu * U[1] * (a * a * a * a * B_term * B_term * B_term * B_term - 1) / (r * r) * f_function * V;      // source term
                // Term 16
                Acc -= nu * U[0] * CoeffUR(r) * f_function * GradV[1];         // first coupling term
                // Term 22 
                Acc -= nu * a * B_term / r * U[2] * (GradV[0] * f_function * CoeffUETA(r) + V * df_function * CoeffUETA(r) + V * f_function * CoeffUETA_deriv(r));    // second couling term
                //Acc -= nu * a * B_term /r U[2] * GradV[0] * f_function * CoeffUETA(r)  
                //Term 26
                //Acc -= nu * a * B_term /r U[2] * V * df_function * CoeffUETA(r) 
                //Term 25 
                //Acc -= nu * a * B_term /r U[2] *  V * f_function * CoeffUETA_deriv(r) 
                //Term 24 
            }

            if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomXI) != 0) {


                double[,] CoeffOfGradUXI = CoeffGradUXI(r);

                for (int d = 0; d < cpv.D; d++) {

                    //Acc -= CoeffOfGradUR[0, d] * GradU[0, d] * (Metric * GradV[0] + V * Metric_dr)
                    //       + CoeffOfGradUR[1, d] * GradU[0, d] * Metric * GradV[1];


                    Acc -= CoeffOfGradUXI[0, d] * GradU[1, d] * GradV[0]
                           + CoeffOfGradUXI[1, d] * GradU[1, d] * GradV[1];
                    // Term 13 & 14

                }
                Acc -= nu * GradU[1, 0] * df_function * V;
                // Term 15

            }
            return (-1.0) * Acc; // Weil NS=Rhs. Terme werden auf die andere Seite der Gleichung gebracht.
        }



        // Inner Edge Form
        public double InnerEdgeForm(ref CommonParams inp,
        double[] U_IN, double[] U_OT, double[,] GradU_IN, double[,] GradU_OT,
        double V_IN, double V_OT, double[] GradV_IN, double[] GradV_OT) {


            double gammaUXI = m_ComputePenalty(m_penaltyFactor, inp.jCellIn, inp.jCellOut, ((GridData)inp.GridDat).Cells.cj);

            double Acc = 0;
            double r = inp.X[0];
            double B_term = Globals.B_term_(r);
            double a = Globals.a;

            double f_function = Globals.f_function_(r);
            double df_function = Globals.df_function_(r);
            double nu = Globals.nu;



            if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomXI) != 0) {
                Acc += nu * 0.5 * (U_IN[1] + U_OT[1]) * f_function / r * inp.Normal[0] * (V_IN - V_OT);
                // Term 17
                Acc += nu * 0.5 * (U_IN[0] + U_OT[0]) * CoeffUR(r) * f_function * inp.Normal[1] * (V_IN - V_OT);
                // first coupling term  
                // Term 21
                Acc += nu * a * B_term / r * 0.5 * (U_IN[2] + U_OT[2]) * CoeffUETA(r) * f_function * inp.Normal[0] * (V_IN - V_OT);   // second coupling term
                // Term 23
                //     
            }


            if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomXI) != 0) {

                double[,] CoeffOfGradUXI = CoeffGradUXI(r);


                // consistency term
                Acc += (CoeffOfGradUXI[0, 0] * 0.5 * (GradU_IN[1, 0] + GradU_OT[1, 0]) * inp.Normal[0]
                    + CoeffOfGradUXI[1, 1] * 0.5 * (GradU_IN[1, 1] + GradU_OT[1, 1]) * inp.Normal[1]) * (V_IN - V_OT);
                // Term 7 & 8



                // symmetry term
                Acc += (CoeffOfGradUXI[0, 0] * 0.5 * (GradV_IN[0] + GradV_OT[0]) * inp.Normal[0]
                       + CoeffOfGradUXI[1, 1] * 0.5 * (GradV_IN[1] + GradV_OT[1]) * inp.Normal[1]) * (U_IN[1] - U_OT[1]);
                // Term 9 & 10
                // penalty term

                Acc -= gammaUXI * (U_IN[1] - U_OT[1]) * (V_IN - V_OT) * Globals.penaltyScaling(r);
                // Wieso kein f(r)? Vorzeichen egal? 
                // Eventuell zählt nur der Absolutwert, daher ist f(r) unwichtig
            }



            return (-1.0) * Acc; // Weil NS=Rhs. Terme werden auf die andere Seite der Gleichung gebracht.
        }



        // Boundary Edge Form

        double[] UDiri(double[] X) {
            double r = X[0];
            double phi = X[1];



            return new double[] {
                Globals.DirichletValue_uR(X),
                Globals.DirichletValue_uXi(X),
                Globals.DirichletValue_uEta(X),
            };
        }
        // now using CoeffGradUR_BE = CoeffGradUR * Metric ,  such that no singularities at r=0 appear
        double[,] CoeffGradUXI_BE(double r) {
            //return new double[,] { { r*r, 0 }, { 0, 1 }};

            double nu = Globals.nu;

            switch (Globals.coordSys) {
                case Globals.CoordSys.pol:
                    if (Globals.AtZeroRadius(r)) {
                        if ((Globals.f_function_(0.1) - 0.01).Abs() > 1.0e-10)
                            throw new ArithmeticException("Someone changed the multiplyer!");

                        return new double[,] { { r * r, 0 }, { 0, 1.0 } };         // multiplier is f(r)=r^2, such that singularity vanishes
                    } else {
                        return new double[,] { { Globals.f_function_(r), 0 }, { 0, Globals.f_function_(r) / (r * r) } };
                    }
                case Globals.CoordSys.hel:
                    if (Globals.AtZeroRadius(r)) {
                        return new double[,] { { nu * Globals.f_function_(r), 0 }, { 0, 1.0 } };         // multiplier is f(r)=B(r)^2, such that singularity vanishes
                    } else {
                        return new double[,] { { nu * Globals.f_function_(r), 0 }, { 0, nu * Globals.f_function_(r) / (Globals.B_term_(r) * Globals.B_term_(r)) } };
                    }
                case Globals.CoordSys.cart:
                    return new double[,] { { 1, 0 }, { 0, 1 } };
                default:
                    throw new NotImplementedException("missing coordinate system " + Globals.coordSys);
            }

        }


        public double BoundaryEdgeForm(ref CommonParamsBnd inp,
            double[] U_IN, double[,] GradU_IN, double V_IN, double[] GradV_IN) {


            double gammaUXI = m_ComputePenalty(m_penaltyFactor, inp.jCellIn, inp.jCellIn, ((GridData)inp.GridDat).Cells.cj);


            double Acc = 0;


            double[] UD;
            UD = UDiri(inp.X);


            double r = inp.X[0];
            double xi = inp.X[1];

            double B_term = Globals.B_term_(r);
            double a = Globals.a;
            double b = Globals.b;
            double nu = Globals.nu;

            if (Globals.AtZeroRadius(r)) {
                //
                // inner edge
                //
                if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomXI) != 0) {

                    //Acc += nu * 0.5 * (U_IN[1] + UD[1]) * Globals.f_function_(r) / r * inp.Normale[0] * V_IN;
                    Acc += nu * U_IN[1] * r / (a * a * r * r + b * b) * inp.Normal[0] * V_IN;
                    // Term 17
                    //Acc += nu * 0.5 * (U_IN[0] + UD[0]) * CoeffUR(r) * Globals.f_function_(r) * inp.Normale[1] * V_IN;      // first coupling term 
                    Acc += nu * U_IN[0] * CoeffURbnd(r) * inp.Normal[1] * V_IN;      // first coupling term
                    // Term 21 
                    //Acc += nu * a * B_term / r * 0.5 * (U_IN[2] + UD[2]) * inp.Normale[0] * CoeffUETA(r) * Globals.f_function_(r) * V_IN;   // second coupling term
                    Acc += nu * a / (Math.Sqrt(a * a * r * r + b * b)) * U_IN[2] * inp.Normal[0] * CoeffUETAbnd(r) * V_IN;   // second coupling term
                    // Term 23
                }

                if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomXI) != 0) {

                    double[,] CoeffOfGradUXI_BE = CoeffGradUXI_BE(r);


                    // consistency term
                    Acc += (CoeffOfGradUXI_BE[0, 0] * (GradU_IN[1, 0]) * inp.Normal[0] + CoeffOfGradUXI_BE[1, 1] * (GradU_IN[1, 1]) * inp.Normal[1]) * (V_IN);
                    // Term 7 & 8

                    // symmetry term
                    // == 0, weil [[u]] = 0;
                    //Acc += (CoeffOfGradUXI_BE[0, 0] * (GradV_IN[0]) * inp.Normale[0] + CoeffOfGradUXI_BE[1, 1] * (GradV_IN[1]) * inp.Normale[1]) * (U_IN[1] - UD[1]);

                    // penalty term
                    // == 0, weil [[u]] = 0;

                    // Ist das die CenterLineCondition?
                    // Wieso ist [[u]] = 0;

                }




            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Dirichlet) {



                if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomXI) != 0) {
                    Acc += nu * 0.5 * (U_IN[1] + UD[1]) * Globals.f_function_(r) / r * inp.Normal[0] * V_IN;
                    // Term 17
                    Acc += nu * 0.5 * (U_IN[0] + UD[0]) * CoeffUR(r) * Globals.f_function_(r) * inp.Normal[1] * V_IN;      // first coupling term  
                                                                                                                           // Term 21
                    Acc += nu * a * B_term / r * 0.5 * (U_IN[2] + UD[2]) * inp.Normal[0] * CoeffUETA(r) * Globals.f_function_(r) * V_IN;   // second coupling term
                    // Term 23
                }
                if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomXI) != 0) {


                    double[,] CoeffOfGradUXI_BE = CoeffGradUXI_BE(r);
                    // consistency term
                    Acc += (CoeffOfGradUXI_BE[0, 0] * (GradU_IN[1, 0]) * inp.Normal[0] +
                        CoeffOfGradUXI_BE[1, 1] * (GradU_IN[1, 1]) * inp.Normal[1]) * (V_IN);
                    // Term 7 & 8

                    // symmetry term
                    Acc += (CoeffOfGradUXI_BE[0, 0] * (GradV_IN[0]) * inp.Normal[0] +
                        CoeffOfGradUXI_BE[1, 1] * (GradV_IN[1]) * inp.Normal[1]) * (U_IN[1] - UD[1]);
                    // Term 9 & 10
                    // Wieso wird der Spung von U genommen, wegen Dirichlet Bedingung?

                    // penalty 
                    Acc -= gammaUXI * (U_IN[1] - UD[1]) * (V_IN) * Globals.penaltyScaling(r);
                    // Term 12

                }

            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Neumann) {
                Acc += 0;
            } else {
                throw new NotImplementedException();
            }
            return (-1.0) * Acc; // Weil NS=Rhs. Terme werden auf die andere Seite der Gleichung gebracht.
        }
    }
}







