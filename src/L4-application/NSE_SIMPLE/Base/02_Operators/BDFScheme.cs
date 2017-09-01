/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;

namespace NSE_SIMPLE {

    /// <summary>
    /// BDF schemes of order 1 to 4.
    /// </summary>
    public class BDFScheme {

        // Coefficients for BDF schemes, cf. Table 1 in 
        // B. Klein, F. Kummer, M. Keil, and M. Oberlack,
        // An extension of the SIMPLE based discontinuous Galerkin solver to unsteady incompressible flows, J. Comput. Phys., 2013.
        public double[] gamma {
            get;
            private set;
        }

        public double[][] beta {
            get;
            private set;
        }

        /// <summary>
        /// Ctor.
        /// </summary>
        public BDFScheme() {

            gamma = new double[4];
            beta = new double[4][];
            beta[0] = new double[2];
            beta[1] = new double[3];
            beta[2] = new double[4];
            beta[3] = new double[5];

            // BDF-1 - Implicit Euler
            gamma[0] = 1.0;
            beta[0][0] = 1.0;
            beta[0][1] = -1.0;

            // BDF-2
            gamma[1] = 2.0;
            beta[1][0] = 3.0;
            beta[1][1] = -4.0;
            beta[1][2] = 1.0;

            // BDF-3
            gamma[2] = 6.0;
            beta[2][0] = 11.0;
            beta[2][1] = -18.0;
            beta[2][2] = 9.0;
            beta[2][3] = -2.0;

            // BDF-4
            gamma[3] = 12.0;
            beta[3][0] = 25.0;
            beta[3][1] = -48.0;
            beta[3][2] = 36.0;
            beta[3][3] = -16.0;
            beta[3][4] = 3.0;
        }

        /// <summary>
        /// Factor for new time step.
        /// </summary>
        /// <returns>
        /// beta_0 / gamma * dt
        /// </returns>
        public double GetLhsSummand(double dt, int BDFOrder) {
            return (beta[BDFOrder - 1][0] / (gamma[BDFOrder - 1] * dt));
        }

        /// <summary>
        /// [Incompressible] Summand for previous time steps in momentum equation.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="BDFOrder"></param>
        /// <param name="Velocity"></param>
        /// <param name="SpatialComponent">Velocity component.</param>
        /// <param name="RhsSummand">Accumulator for the result.</param>
        public void ComputeRhsSummand(double dt, int BDFOrder,
            VectorFieldHistory<SinglePhaseField> Velocity, int SpatialComponent,
            SinglePhaseField RhsSummand) {

            for (int alpha = 1; alpha <= BDFOrder; alpha++)
                RhsSummand.Acc(beta[BDFOrder - 1][alpha], Velocity[1 - alpha][SpatialComponent]);

            RhsSummand.Scale(1.0 / (gamma[BDFOrder - 1] * dt));
        }

        /// <summary>
        /// [LowMach] Summand for previous time steps in momentum equation.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="BDFOrder"></param>
        /// <param name="Scalar"></param>
        /// <param name="Velocity"></param>
        /// <param name="SpatialComponent">Velocity component.</param>
        /// <param name="EoS"></param>
        /// <param name="RhsSummand">Accumulator for the result.</param>
        public void ComputeRhsSummand(double dt, int BDFOrder,
            ScalarFieldHistory<SinglePhaseField> Scalar, VectorFieldHistory<SinglePhaseField> Velocity, int SpatialComponent, MaterialLaw EoS,
            SinglePhaseField RhsSummand) {

            for (int alpha = 1; alpha <= BDFOrder; alpha++)
                RhsSummand.ProjectFunction(beta[BDFOrder - 1][alpha],
                    (X, U, cell) => EoS.GetDensity(U[0]) * U[1],
                    null,
                    Scalar[1 - alpha],
                    Velocity[1 - alpha][SpatialComponent]);

            RhsSummand.Scale(1.0 / (gamma[BDFOrder - 1] * dt));
        }

        /// <summary>
        /// [LowMach] Summand for previous time steps in temperature equation.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="BDFOrder"></param>
        /// <param name="Temperature"></param>
        /// <param name="EoS"></param>
        /// <param name="RhsSummand">Accumulator for the result.</param>
        public void ComputeRhsSummand(double dt, int BDFOrder,
            ScalarFieldHistory<SinglePhaseField> Temperature, MaterialLaw EoS,
            SinglePhaseField RhsSummand) {

            for (int alpha = 1; alpha <= BDFOrder; alpha++)
                RhsSummand.ProjectFunction(beta[BDFOrder - 1][alpha],
                    (X, U, cell) => EoS.GetDensity(U[0]) * U[0],
                    null,
                    Temperature[1 - alpha]);

            RhsSummand.Scale(1.0 / (gamma[BDFOrder - 1] * dt));
        }

        /// <summary>
        /// [Multiphase] Summand for previous time steps in level-set equation.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="BDFOrder"></param>        
        /// <param name="Scalar"></param>        
        /// <param name="RhsSummand">Accumulator for the result.</param>
        public void ComputeRhsSummand(double dt, int BDFOrder,
            ScalarFieldHistory<SinglePhaseField> Scalar,
            SinglePhaseField RhsSummand) {

            for (int alpha = 1; alpha <= BDFOrder; alpha++)
                RhsSummand.Acc(beta[BDFOrder - 1][alpha], Scalar[1 - alpha]);

            RhsSummand.Scale(1.0 / (gamma[BDFOrder - 1] * dt));
        }

        /// <summary>
        /// [LowMach] Summand of all time steps for scalar variables, which are constant in space.
        /// Used for time derivative of thermodynamic pressure in Low-Mach flows.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="BDFOrder"></param>
        /// <param name="Scalar"></param>
        /// <returns>
        /// Summand of all time steps of <paramref name="Scalar"/>.
        /// </returns>
        //public double ComputeSummandScalarHistory(double dt, int BDFOrder, ScalarFieldHistory<SinglePhaseField> Scalar) {
        //    double Summand = 0.0;

        //    for (int alpha = 0; alpha <= BDFOrder; alpha++)
        //        Summand += beta[BDFOrder - 1][alpha] * Scalar[1 - alpha].GetMeanValue(0);

        //    Summand *= 1.0 / (gamma[BDFOrder - 1] * dt);

        //    return Summand;
        //}

        /// <summary>
        /// [LowMach] Summand of all time steps for density.
        /// Used for time derivative of density in Corrector for Low-Mach flows.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="BDFOrder"></param>
        /// <param name="Temperature"></param>
        /// <param name="EoS"></param>
        /// <param name="RhsSummand"></param>
        public void ComputeDensitySummand(double dt, int BDFOrder,
            ScalarFieldHistory<SinglePhaseField> Temperature, MaterialLaw EoS,
            SinglePhaseField RhsSummand) {

            for (int alpha = 0; alpha <= BDFOrder; alpha++)
                RhsSummand.ProjectFunction(beta[BDFOrder - 1][alpha],
                    (X, U, cell) => EoS.GetDensity(U[0]),
                    null,
                    Temperature[1 - alpha]);

            RhsSummand.Scale(1.0 / (gamma[BDFOrder - 1] * dt));
        }
    }
}