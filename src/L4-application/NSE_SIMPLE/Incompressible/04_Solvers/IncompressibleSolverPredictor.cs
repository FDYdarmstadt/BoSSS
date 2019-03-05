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
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using BoSSS.Solution;

namespace NSE_SIMPLE {
    /// <summary>
    /// Predictor solver for incompressible flows.
    /// </summary>
    public class SolverPredictor : SIMPLESolver {

        SolverConfiguration m_SolverConf;

        SIMPLEMatrixAssembly[] m_MatAsmblyPredictor;
        SIMPLEMatrixAssembly m_MatAsmblyPredictorApprox;
        SIMPLEOperator[] m_PressureGradient;

        // Time discretization
        BDFScheme m_BDF = null;

        VectorFieldHistory<SinglePhaseField> m_Velocity;
        SinglePhaseField m_Pressure;

        double m_RelaxFactor;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="solverConf"></param>
        /// <param name="_sparseSolver"></param>
        /// <param name="MatAsmblyPredictor"></param>
        /// <param name="MatAsmblyPredictorApprox"></param>
        /// <param name="PressureGradient"></param>
        /// <param name="BDF"></param>
        /// <param name="Velocity"></param>
        /// <param name="Pressure"></param>
        public SolverPredictor(SolverConfiguration solverConf, ISparseSolver _sparseSolver,
            SIMPLEMatrixAssembly[] MatAsmblyPredictor, SIMPLEMatrixAssembly MatAsmblyPredictorApprox, SIMPLEOperator[] PressureGradient, BDFScheme BDF,
            VectorFieldHistory<SinglePhaseField> Velocity, SinglePhaseField Pressure)
            : base(solverConf, _sparseSolver) {

            m_SolverConf = solverConf;

            m_MatAsmblyPredictor = MatAsmblyPredictor;
            m_MatAsmblyPredictorApprox = MatAsmblyPredictorApprox;
            m_PressureGradient = PressureGradient;

            m_BDF = BDF;

            m_Velocity = Velocity;
            m_Pressure = Pressure;

            m_RelaxFactor = (1.0 - base.m_solverConf.Control.RelexationFactorVelocity) / base.m_solverConf.Control.RelexationFactorVelocity;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="dt"></param>
        /// <returns></returns>
        protected override MsrMatrix DefineMatrix(double dt) {
            MsrMatrix res = new MsrMatrix(m_MatAsmblyPredictor[0].AssemblyMatrix);

            //See left-hand side of Eq. (18) in 
            //B. Klein, F. Kummer, M. Keil, and M. Oberlack,
            //An extension of the SIMPLE based discontinuous Galerkin solver to unsteady incompressible flows, J. Comput. Phys., 2013.
            if (m_RelaxFactor != 0.0)
                res.Acc(m_RelaxFactor, m_MatAsmblyPredictorApprox.AssemblyMatrix);

            if (m_BDF != null) {
                double LhsSummand = m_BDF.GetLhsSummand(dt, base.m_solverConf.BDFOrder);
                res.AccEyeSp(LhsSummand);
            }

            return res;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="SpatialComponent"></param>
        /// <returns></returns>
        protected override IList<double> DefineRhs(double dt, int SpatialComponent) {
            double[] rhs = new double[m_MatAsmblyPredictor[SpatialComponent].LocalLength];

            double[] PredictorAffine = m_MatAsmblyPredictor[SpatialComponent].AssemblyAffine;

            SinglePhaseField RhsSummand = new SinglePhaseField(m_Velocity.Current[0].Basis);
            if (m_BDF != null)
                m_BDF.ComputeRhsSummand(dt, base.m_solverConf.BDFOrder, m_Velocity, SpatialComponent, RhsSummand);

            double[] PressureGradient = new double[rhs.Length];
            m_PressureGradient[SpatialComponent].OperatorMatrix.SpMVpara(1.0, m_Pressure.CoordinateVector, 0.0, PressureGradient);

            double[] UnderRelaxation = new double[rhs.Length];
            if (m_RelaxFactor != 0.0)
                m_MatAsmblyPredictorApprox.AssemblyMatrix.SpMVpara(m_RelaxFactor,
                    m_Velocity.Current[SpatialComponent].CoordinateVector,
                    0.0,
                    UnderRelaxation);

            for (int i = 0; i < rhs.Length; i++) {
                //See right-hand side of Eq. (18) in 
                //B. Klein, F. Kummer, M. Keil, and M. Oberlack,
                //An extension of the SIMPLE based discontinuous Galerkin solver to unsteady incompressible flows, J. Comput. Phys., 2013.
                rhs[i] =
                    -PredictorAffine[i]
                    - m_PressureGradient[SpatialComponent].OperatorAffine[i]
                    - RhsSummand.CoordinateVector[i]
                    - PressureGradient[i]
                    + UnderRelaxation[i];
            }

            return rhs;
        }
    }
}
