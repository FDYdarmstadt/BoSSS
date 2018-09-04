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
using ilPSP.Utils;
using BoSSS.Foundation;
using BoSSS.Solution;
using BoSSS.Platform;
using BoSSS.Solution.NSECommon;

namespace NSE_SIMPLE {

    /// <summary>
    /// Predictor solver for variable density flows.
    /// </summary>
    public class VariableDensitySolverPredictor : SIMPLESolver {

        BlockDiagonalMatrix m_DensityMatrix;
        SIMPLEMatrixAssembly m_MatAsmblyPredictor;
        SIMPLEMatrixAssembly m_MatAsmblyPredictorApprox;

        SIMPLEOperator[] m_PressureGradient;
        SinglePhaseField m_Pressure;
        SIMPLEMatrixAssembly[,] m_MatAsmblyViscSplit;
        IEvaluatorNonLin[] m_BuoyantForceEvaluator;
        VectorFieldHistory<SinglePhaseField> m_Velocity;
        ScalarFieldHistory<SinglePhaseField> m_Scalar;
        MaterialLaw m_EoS;

        double m_RelaxFactor;

        // Time discretization
        BDFScheme m_BDF = null;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="solverConf"></param>
        /// <param name="sparseSolver"></param>
        /// <param name="DensityMatrix"></param>
        /// <param name="MatAsmblyPredictor"></param>
        /// <param name="MatAsmblyPredictorApprox"></param>
        /// <param name="PressureGradient"></param>
        /// <param name="Pressure"></param>
        /// <param name="MatAsmblyViscSplit"></param>
        /// <param name="BuoyantForce"></param>
        /// <param name="Velocity"></param>
        /// <param name="Scalar"></param>
        /// <param name="EoS"></param>
        /// <param name="BDF"></param>
        public VariableDensitySolverPredictor(SolverConfiguration solverConf, ISparseSolver sparseSolver,
            BlockDiagonalMatrix DensityMatrix, SIMPLEMatrixAssembly MatAsmblyPredictor, SIMPLEMatrixAssembly MatAsmblyPredictorApprox,
            SIMPLEOperator[] PressureGradient, SinglePhaseField Pressure,
            SIMPLEMatrixAssembly[,] MatAsmblyViscSplit, IEvaluatorNonLin[] BuoyantForce,
            VectorFieldHistory<SinglePhaseField> Velocity, ScalarFieldHistory<SinglePhaseField> Scalar, MaterialLaw EoS, BDFScheme BDF)
            : base(solverConf, sparseSolver) {

            m_DensityMatrix = DensityMatrix;
            m_MatAsmblyPredictor = MatAsmblyPredictor;
            m_MatAsmblyPredictorApprox = MatAsmblyPredictorApprox;

            m_PressureGradient = PressureGradient;
            m_Pressure = Pressure;
            m_MatAsmblyViscSplit = MatAsmblyViscSplit;
            m_BuoyantForceEvaluator = BuoyantForce;
            m_Velocity = Velocity;
            m_Scalar = Scalar;
            m_EoS = EoS;

            m_RelaxFactor = (1.0 - base.m_solverConf.Control.RelexationFactorVelocity) / base.m_solverConf.Control.RelexationFactorVelocity;

            m_BDF = BDF;
        }

        protected override MsrMatrix DefineMatrix(double dt) {
            MsrMatrix res = new MsrMatrix(m_MatAsmblyPredictor.AssemblyMatrix);

            if (m_RelaxFactor != 0.0)
                res.Acc(m_RelaxFactor, m_MatAsmblyPredictorApprox.AssemblyMatrix);

            if (m_BDF != null) {
                double LhsSummand = m_BDF.GetLhsSummand(dt, base.m_solverConf.BDFOrder);
                res.Acc(LhsSummand, m_DensityMatrix);
            }

            return res;
        }

        protected override IList<double> DefineRhs(double dt, int SpatialComponent) {
            double[] rhs = new double[m_MatAsmblyPredictor.LocalLength];

            double[] PredictorAffine = m_MatAsmblyPredictor.AssemblyAffine;

            SinglePhaseField RhsSummand = new SinglePhaseField(m_Velocity.Current[SpatialComponent].Basis);
            if (m_BDF != null)
                m_BDF.ComputeRhsSummand(dt, base.m_solverConf.BDFOrder, m_Scalar, m_Velocity, SpatialComponent, m_EoS, RhsSummand);

            double[] PressureGradient = new double[rhs.Length];
            m_PressureGradient[SpatialComponent].OperatorMatrix.SpMVpara(1.0, m_Pressure.CoordinateVector, 0.0, PressureGradient);

            double[] ViscSplitExplicitPart = new double[rhs.Length];
            for (int j = 0; j < m_Velocity.Current.Dim; j++) {
                if (j != SpatialComponent) {
                    m_MatAsmblyViscSplit[SpatialComponent, j].AssemblyMatrix.SpMVpara(1.0, m_Velocity.Current[j].CoordinateVector, 1.0, ViscSplitExplicitPart);
                }
            }

            double[] BuoyantForce = new double[rhs.Length];
            if (m_BuoyantForceEvaluator != null) {
                m_BuoyantForceEvaluator[SpatialComponent].Evaluate(1.0, 0.0, BuoyantForce);
            }

            double[] UnderRelaxation = new double[rhs.Length];
            if (m_RelaxFactor != 0.0)
                m_MatAsmblyPredictorApprox.AssemblyMatrix.SpMVpara(m_RelaxFactor, m_Velocity.Current[SpatialComponent].CoordinateVector, 0.0, UnderRelaxation);


            for (int i = 0; i < rhs.Length; i++) {
                rhs[i] =
                    -PredictorAffine[i]
                    - m_PressureGradient[SpatialComponent].OperatorAffine[i]
                    - RhsSummand.CoordinateVector[i]
                    - PressureGradient[i]
                    - ViscSplitExplicitPart[i]
                    + BuoyantForce[i]
                    + UnderRelaxation[i];
            }

            return rhs;
        }
    }
}
