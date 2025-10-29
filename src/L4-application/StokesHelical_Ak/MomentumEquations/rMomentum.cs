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

namespace StokesHelical_Ak {
    class rMomentum :
            IEdgeForm, // edge integrals
            IVolumeForm     // volume integrals

    {



        public TermSwitch m_TermSwitch;



        private Func<double, int, int, MultidimensionalArray, double> m_ComputePenalty;

        private double m_penaltyFactor;


        public rMomentum(TermSwitch termSwitch, double PenaltyFactor, Func<double, int, int, MultidimensionalArray, double> ComputePenalty = null) {
            m_TermSwitch = termSwitch;
            m_ComputePenalty = ComputePenalty;
            m_penaltyFactor = PenaltyFactor;
        }




        /// The parameter list for the divergence is empty:
        public IList<string> ParameterOrdering {
            get { return null; }
        }


        // Velocity components

        public IList<String> ArgumentOrdering {
            get { return new string[] { "ur", "uxi", "ueta" }; }
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


        double[,] CoeffGradUR(double r) {

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

            double a = Globals.a;
            double b = Globals.b;
            double B_term = Globals.B_term_(r);
            double nu = Globals.nu;





            if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomR) != 0) {
                Acc -= nu * U[0] * ((df_function * r - f_function) / (r * r) * V + f_function / r * GradV[0]);
                // Term 17, Acc -= nu *   U[0] * V * df_funciton r / (r * r)
                // Term 18, Acc -= nu * - U[0] * V *  f_function / (r * r)
                // Term 19, Acc -= nu *   U[0] * GradV[0] * f_function / r 
                Acc -= nu * 1 / (r * r) * U[0] * f_function * V;   // source term in diffusion operator
                // Term 15
                Acc -= nu * U[2] * (-2 * a * b * B_term / (r * r) * f_function) * GradV[1];       // coupling term ueta
                // Term 21
                Acc -= nu * U[1] * (-2 * b * b * B_term / (r * r * r) * f_function) * GradV[1];   // coupling term uxi
                // Term 23
            }

            if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomR) != 0) {


                double[,] CoeffOfGradUR = CoeffGradUR(r);

                for (int d = 0; d < cpv.D; d++) {


                    Acc -= CoeffOfGradUR[0, d] * GradU[0, d] * GradV[0]
                           + CoeffOfGradUR[1, d] * GradU[0, d] * GradV[1];  // Term 12 und 13


                }
                Acc -= nu * GradU[0, 0] * df_function * V; // Term 14

            }
            return (-1.0) * Acc; // Multiplikation mit -1, da die NS-Gleichung=RHS
        }



        // Inner Edge Form
        public double InnerEdgeForm(ref CommonParams inp,
        double[] U_IN, double[] U_OT, double[,] GradU_IN, double[,] GradU_OT,
        double V_IN, double V_OT, double[] GradV_IN, double[] GradV_OT) {


            double gammaUR = m_ComputePenalty(m_penaltyFactor, inp.jCellIn, inp.jCellOut, ((GridData)inp.GridDat).Cells.cj);


            double Acc = 0;
            double r = inp.X[0];


            double f_function = Globals.f_function_(r);
            double df_function = Globals.df_function_(r);

            double a = Globals.a;
            double b = Globals.b;
            double B_term = Globals.B_term_(r);
            double nu = Globals.nu;


            if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomR) != 0) {
                Acc += nu * 0.5 * (U_IN[0] + U_OT[0]) * f_function / r * inp.Normal[0] * (V_IN - V_OT);
                // Term 16
                Acc += nu * 0.5 * (U_IN[2] + U_OT[2]) * (-2 * a * b * B_term / (r * r)) * f_function * inp.Normal[1] * (V_IN - V_OT);    // coupling term ueta
                // Term 20
                Acc += nu * 0.5 * (U_IN[1] + U_OT[1]) * (-2 * b * b * B_term / (r * r * r)) * f_function * inp.Normal[1] * (V_IN - V_OT);    // coupling term ueta
                // Term 22
            }

            if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomR) != 0) {

                double[,] CoeffOfGradUR = CoeffGradUR(r);


                // consistency term
                Acc += (CoeffOfGradUR[0, 0] * 0.5 * (GradU_IN[0, 0] + GradU_OT[0, 0]) * inp.Normal[0]
                    + CoeffOfGradUR[1, 1] * 0.5 * (GradU_IN[0, 1] + GradU_OT[0, 1]) * inp.Normal[1]) * (V_IN - V_OT);
                // Term 6 & 7

                // symmetry term
                Acc += (CoeffOfGradUR[0, 0] * 0.5 * (GradV_IN[0] + GradV_OT[0]) * inp.Normal[0]
                       + CoeffOfGradUR[1, 1] * 0.5 * (GradV_IN[1] + GradV_OT[1]) * inp.Normal[1]) * (U_IN[0] - U_OT[0]);
                // Term 8 & 9

                // penalty term
                Acc -= gammaUR * (U_IN[0] - U_OT[0]) * (V_IN - V_OT) * Globals.penaltyScaling(r);
                // Wo ist  f_function ???
                // Term 11
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
        // Also: f(r) = B(r)^2
        double[,] CoeffGradUR_BE(double r) {
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
                        return new double[,] { { nu * Globals.f_function_(r), 0 }, { 0, nu * 1.0 } };         // multiplier is f(r)=B(r)^2, such that singularity vanishes
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


            double gammaUR = m_ComputePenalty(m_penaltyFactor, inp.jCellIn, inp.jCellIn, ((GridData)inp.GridDat).Cells.cj);


            double Acc = 0;


            double[] UD;
            UD = UDiri(inp.X);




            double r = inp.X[0];
            double xi = inp.X[1];

            double f_function = Globals.f_function_(r);

            double a = Globals.a;
            double b = Globals.b;
            double B_term = Globals.B_term_(r);
            double nu = Globals.nu;


            if (Globals.AtZeroRadius(r)) {
                //
                // inner edge
                //
                if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomR) != 0) {

                    //Acc += nu  * U_IN[0]  * f_function / r * inp.Normale[0] * V_IN;
                    Acc += nu * U_IN[0] * r / (a * a * r * r + b * b) * inp.Normal[0] * V_IN;         // multiplicator f(r)=B^2
                                                                                                      // Term 16

                    //Acc += nu * U_IN[2]  * (-2 * a * b * B_term / (r * r)) * f_function * inp.Normale[1] * V_IN;
                    Acc += nu * U_IN[2] * (-2 * a * b * B_term / (a * a * r * r + b * b)) * inp.Normal[1] * V_IN;  // multiplicator f(r)=B^2
                                                                                                                   // Term 20


                    //Acc += nu *  U_IN[1]  * (-2 * b * b * B_term / (r * r * r)) * f_function * inp.Normale[1] * V_IN;
                    Acc += nu * U_IN[1] * (-2 * b * b / ((a * a * r * r + b * b) * Math.Sqrt(a * a * r * r + b * b))) * inp.Normal[1] * V_IN;   // multiplicator f(r)=B^2
                                                                                                                                                // Term 22
                }
                // Nur die Innenwerte, weil Boundary

                if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomR) != 0) {

                    double[,] CoeffOfGradUR_BE = CoeffGradUR_BE(r);

                    // consistency term
                    Acc += (CoeffOfGradUR_BE[0, 0] * (GradU_IN[0, 0]) * inp.Normal[0] + CoeffOfGradUR_BE[1, 1] * (GradU_IN[0, 1]) * inp.Normal[1]) * (V_IN);
                    // Term 6 und 7

                    // symmetry term
                    // == 0, weil [[u]] = 0;
                    // Acc += (CoeffOfGradUR_BE[0, 0] * (GradV_IN[0]) * inp.Normale[0] + CoeffOfGradUR_BE[1, 1] * (GradV_IN[1]) * inp.Normale[1]) * (U_IN[0] - UD[0]);
                    // Wieso ist der Sprung am Innteren Boundary = 0?


                    // penalty term
                    // == 0, weil [[u]] = 0;

                }

            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Dirichlet) {

                // mit Dirichlet Boundary
                if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomR) != 0) {
                    Acc += nu * 0.5 * (U_IN[0] + UD[0]) * Globals.f_function_(r) / r * inp.Normal[0] * V_IN;
                    //Acc += nu * 0.5 * (U_IN[0]+UD[0]) * f_function      / r * inp.Normal[0] * V_IN;
                    //Term 16
                    Acc += nu * 0.5 * (U_IN[2] + UD[2]) * (-2 * a * b * B_term / (r * r)) * f_function * inp.Normal[1] * V_IN;
                    //  Term 20
                    Acc += nu * 0.5 * (U_IN[1] + UD[1]) * (-2 * b * b * B_term / (r * r * r)) * f_function * inp.Normal[1] * V_IN;
                    //  Term 22 
                }

                if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomR) != 0) {

                    double[,] CoeffOfGradUR_BE = CoeffGradUR_BE(r);

                    // consistency term
                    Acc += (CoeffOfGradUR_BE[0, 0] * (GradU_IN[0, 0]) * inp.Normal[0] + CoeffOfGradUR_BE[1, 1] * (GradU_IN[0, 1]) * inp.Normal[1]) * (V_IN);
                    // Term 6 und 7

                    // symmetry term
                    Acc += (CoeffOfGradUR_BE[0, 0] * (GradV_IN[0]) * inp.Normal[0] + CoeffOfGradUR_BE[1, 1] * (GradV_IN[1]) * inp.Normal[1]) * (U_IN[0] - UD[0]);
                    // Term 8 und 9
                    // Wieso Differenz?

                    // penalty term
                    Acc -= gammaUR * (U_IN[0] - UD[0]) * (V_IN) * Globals.penaltyScaling(r);
                    // Wo ist f(r)?
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







