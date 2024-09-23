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
    class etaMomentum :
            IEdgeForm, // edge integrals
            IVolumeForm,     // volume integrals
        ISupportsJacobianComponent // For Jacobian, required for Newton Solver
        {



        public TermSwitch m_TermSwitch;



        private Func<double, int, int, MultidimensionalArray, double> m_ComputePenalty;

        private double m_penaltyFactor;

        int ResolutionXI;


        public etaMomentum(int ResolutionXI, TermSwitch termSwitch, double PenaltyFactor, Func<double, int, int, MultidimensionalArray, double> ComputePenalty = null) {
            m_TermSwitch = termSwitch;
            m_ComputePenalty = ComputePenalty;
            m_penaltyFactor = PenaltyFactor;
            this.ResolutionXI = ResolutionXI;
        }

        /// <summary>
        /// coefficient of durdxi
        /// </summary>
        /// <param name="r"></param>
        /// <returns></returns>
        static double CoeffUR(double r) {
            return 2 * Globals.a * Globals.b * Globals.B_term_(r) / (r * r);
        }

        /// <summary>
        /// coefficient of durdxi at the boundary r=0: this is coefficient of durduxi * f_function(r)
        /// </summary>
        /// <param name="r"></param>
        /// <returns></returns>
        static double CoeffURbnd(double r) {
            double a = Globals.a;
            double b = Globals.b;
            return 2 * Globals.a * Globals.b * Globals.B_term_(r) / (a * a * r * r + b * b);
        }


        /// <summary>
        /// coefficient of duxidr
        /// </summary>
        /// <param name="r"></param>
        /// <returns></returns>
        static double CoeffUXI(double r) {
            return -2.0 * Globals.a * Globals.b * Globals.B_term_(r) / (r * r);
        }

        /// <summary>
        /// coefficient of duxidr at the boundary r=0: this is coefficient of durduxi * f_function(r)
        /// </summary>
        /// <param name="r"></param>
        /// <returns></returns>
        static double CoeffUXIbnd(double r) {
            double a = Globals.a;
            double b = Globals.b;
            return -2.0 * Globals.a * Globals.b * Globals.B_term_(r) / (a * a * r * r + b * b);
        }



        static double CoeffUXI_deriv(double r) {
            double B_term = Globals.B_term_(r);
            double a = Globals.a;
            double b = Globals.b;
            return (2.0 * a * b * B_term / (r * r * r) + 2.0 * a * a * a * b * B_term * B_term * B_term / (r * r * r));
            // 100 % Confirmed!
        }

        /// The parameter list for the divergence is empty:
        public IList<string> ParameterOrdering {
            get { return null; }
        }


        // Velocity components

        public IList<String> ArgumentOrdering {
            get { return new string[] {  "Velocity_R", "Velocity_XI", "Velocity_ETA" , "Pressure" }; }
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


        double[,] CoeffGradUETA(double r) {

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
            if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomETA) != 0) {
                Acc -= nu * U[2] * ((df_function * r - f_function) / (r * r) * V + f_function / r * GradV[0]);
                // Acc -= nu * U[2] * df_function * V * r / (r * r); Term 17
                // Acc -= - nu * U[2] * f_function * V / (r * r); Term 18
                // Acc -= nu * U[2] * f_function / r * GradV[0]; Term 19
                //
                Acc += nu * U[2] * a * a * B_term * B_term * (a * a * B_term * B_term - 2) / (r * r) * f_function * V; // source term
                // Term 15

                Acc -= nu * U[0] * CoeffUR(r) * f_function * GradV[1];     // first coupling term
                // Term 21
                Acc -= nu * B_term * U[1] * (GradV[0] * f_function * CoeffUXI(r) + V * df_function * CoeffUXI(r) + V * f_function * CoeffUXI_deriv(r));     // second coupling term
                // Acc -= nu * B_term * U[1] * GradV[0] * f_function * CoeffUXI(r); Term 25
                // Acc -= nu * B_term * U[1] * V * df_function * CoeffUXI(r); Term 24
                // Acc -= nu * B_term * U[1] * V * f_function * CoeffUXI_deriv(r); Term 23
            }

            if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomETA) != 0) {


                double[,] CoeffOfGradUETA = CoeffGradUETA(r);

                for (int d = 0; d < cpv.D; d++) {


                    Acc -= CoeffOfGradUETA[0, d] * GradU[2, d] * GradV[0]
                           + CoeffOfGradUETA[1, d] * GradU[2, d] * GradV[1];
                    // Term 12 & 13

                }



                // Term from derivative of multiplier f(r)
                Acc -= nu * GradU[2, 0] * df_function * V;
                // Term 14
            }
            return (-1.0) * Acc;
        }



        // Inner Edge Form
        public double InnerEdgeForm(ref CommonParams inp,
        double[] U_IN, double[] U_OT, double[,] GradU_IN, double[,] GradU_OT,
        double V_IN, double V_OT, double[] GradV_IN, double[] GradV_OT) {


            double gammaUETA = m_ComputePenalty(m_penaltyFactor, inp.jCellIn, inp.jCellOut, ((GridData)inp.GridDat).Cells.cj);


            double Acc = 0;
            double r = inp.X[0];

            double B_term = Globals.B_term_(r);
            double f_function = Globals.f_function_(r);
            double df_function = Globals.df_function_(r);
            double nu = Globals.nu;




            if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomETA) != 0) {
                Acc += nu * 0.5 * (U_IN[2] + U_OT[2]) * f_function / r * inp.Normal[0] * (V_IN - V_OT);
                // Term 16
                Acc += nu * 0.5 * (U_IN[0] + U_OT[0]) * CoeffUR(r) * f_function * inp.Normal[1] * (V_IN - V_OT);      // first coupling term  
                // Term 20
                Acc += nu * B_term * 0.5 * (U_IN[1] + U_OT[1]) * CoeffUXI(r) * f_function * inp.Normal[0] * (V_IN - V_OT);   // second coupling term
                // Term 22
            }


            if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomETA) != 0) {

                double[,] CoeffOfGradUETA = CoeffGradUETA(r);
                // consistency term
                Acc += (CoeffOfGradUETA[0, 0] * 0.5 * (GradU_IN[2, 0] + GradU_OT[2, 0]) * inp.Normal[0]
                    + CoeffOfGradUETA[1, 1] * 0.5 * (GradU_IN[2, 1] + GradU_OT[2, 1]) * inp.Normal[1]) * (V_IN - V_OT);
                // Term 6 & 7


                // symmetry term
                Acc += (CoeffOfGradUETA[0, 0] * 0.5 * (GradV_IN[0] + GradV_OT[0]) * inp.Normal[0]
                       + CoeffOfGradUETA[1, 1] * 0.5 * (GradV_IN[1] + GradV_OT[1]) * inp.Normal[1]) * (U_IN[2] - U_OT[2]);
                // Term 8 & 9

                // penalty term

                //Acc -= gammaUETA * (U_IN[2] - U_OT[2]) * (V_IN - V_OT) * Math.Max(1.0, Math.Abs(Globals.f_function_(r) * Globals.f_function_(r) / (r * r)));
                Acc -= gammaUETA * (U_IN[2] - U_OT[2]) * (V_IN - V_OT) * Globals.penaltyScaling(r);
                // Term 11

            }
            // pressure stabilization
            if (Globals.pressureStabilEtaMom == true) {
                Acc -= 1 / Globals.nu * 2 * Math.PI / ResolutionXI * (U_IN[3] - U_OT[3]) * (V_IN - V_OT) * Globals.penaltyScaling(r);
                // All Dead Code! Never true. See Details on DiPietro 2013
            }
            return (-1.0) * Acc;
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
        double[,] CoeffGradUETA_BE(double r) {
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


            double gammaUETA = m_ComputePenalty(m_penaltyFactor, inp.jCellIn, inp.jCellIn, ((GridData)inp.GridDat).Cells.cj);


            double Acc = 0;


            double[] UD;
            UD = UDiri(inp.X);

            double r = inp.X[0];
            double B_term = Globals.B_term_(r);
            double xi = inp.X[1];


            double nu = Globals.nu;
            double a = Globals.a;
            double b = Globals.b;

            if (Globals.AtZeroRadius(r)) {
                //
                // inner edge
                //
                if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomETA) != 0) {


                    //Acc += nu * 0.5 * (U_IN[2] + UD[2]) * Globals.f_function_(r) / r * inp.Normale[0] * V_IN;
                    Acc += nu * U_IN[2] * r / (a * a * r * r + b * b) * inp.Normal[0] * V_IN;
                    // Term 16

                    //Acc += nu * 0.5 * (U_IN[0] + UD[0]) * CoeffUR(r) * Globals.f_function_(r) * inp.Normale[1] * V_IN;      // first coupling term  
                    Acc += nu * U_IN[0] * CoeffURbnd(r) * inp.Normal[1] * V_IN;      // first coupling term  
                    // Term 20

                    //Acc += nu * B_term * 0.5 * (U_IN[1] + UD[1]) * inp.Normale[0] * CoeffUXI(r) * Globals.f_function_(r) * V_IN;   // second coupling term
                    Acc += nu * B_term * U_IN[1] * inp.Normal[0] * CoeffUXIbnd(r) * V_IN;   // second coupling term
                    // Term 22
                }
                if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomETA) != 0) {

                    double[,] CoeffOfGradUETA_BE = CoeffGradUETA_BE(r);


                    // consistency term
                    Acc += (CoeffOfGradUETA_BE[0, 0] * (GradU_IN[2, 0]) * inp.Normal[0] + CoeffOfGradUETA_BE[1, 1] * (GradU_IN[2, 1]) * inp.Normal[1]) * (V_IN);
                    // Term 6 & 7

                    // symmetry term
                    // == 0, weil [[u]] = 0;
                    //Acc += (CoeffOfGradUETA_BE[0, 0] * (GradV_IN[0]) * inp.Normale[0] + CoeffOfGradUETA_BE[1, 1] * (GradV_IN[1]) * inp.Normale[1]) * (U_IN[2] - UD[2]);



                    // penalty term
                    // == 0, weil [[u]] = 0;

                }


            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Dirichlet) {


                if ((this.m_TermSwitch & TermSwitch.Viscosity_1stOrder_MomETA) != 0) {

                    Acc += nu * 0.5 * (U_IN[2] + UD[2]) * Globals.f_function_(r) / r * inp.Normal[0] * V_IN;
                    // Term 16
                    Acc += nu * 0.5 * (U_IN[0] + UD[0]) * CoeffUR(r) * Globals.f_function_(r) * inp.Normal[1] * V_IN;      // first coupling term
                    // Term 20                                                                                              
                    Acc += nu * B_term * 0.5 * (U_IN[1] + UD[1]) * inp.Normal[0] * CoeffUXI(r) * Globals.f_function_(r) * V_IN;   // second coupling term
                    // Term 22

                }
                if ((this.m_TermSwitch & TermSwitch.Viscosity_2ndOrder_MomETA) != 0) {

                    double[,] CoeffOfGradUETA_BE = CoeffGradUETA_BE(r);


                    // consistency term
                    Acc += (CoeffOfGradUETA_BE[0, 0] * (GradU_IN[2, 0]) * inp.Normal[0] + CoeffOfGradUETA_BE[1, 1] * (GradU_IN[2, 1]) * inp.Normal[1]) * (V_IN);
                    // Term 6 & 7

                    // symmetry term
                    Acc += (CoeffOfGradUETA_BE[0, 0] * (GradV_IN[0]) * inp.Normal[0] + CoeffOfGradUETA_BE[1, 1] * (GradV_IN[1]) * inp.Normal[1]) * (U_IN[2] - UD[2]);
                    // Term 8 & 9
                    // penalty term
                    //Acc -= gammaUETA * (U_IN[2] - UD[2]) * (V_IN) * Math.Max(1.0, Math.Abs(Globals.f_function_(r) * Globals.f_function_(r) / (r * r)));
                    Acc -= gammaUETA * (U_IN[2] - UD[2]) * (V_IN) * Globals.penaltyScaling(r);
                    // Term 11
                }
                // pressure stabilization
                if (Globals.pressureStabilEtaMom == true) {
                    Acc -= 1 / Globals.nu * 2 * Math.PI / ResolutionXI * (U_IN[3] - UD[3]) * (V_IN) * Globals.penaltyScaling(r);
                    // All Dead Code! Never true. See Details on DePietro 2013
                }


            } else if (Globals.BoundaryType(inp.X) == BoundaryTypeE.Neumann) {
                Acc += 0;
            } else {
                throw new NotImplementedException();
            }
            return (-1.0) * Acc;
        }

        /// <summary>
        /// Linear component - derivative is just this.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new[] { this };
        }
    }
}







