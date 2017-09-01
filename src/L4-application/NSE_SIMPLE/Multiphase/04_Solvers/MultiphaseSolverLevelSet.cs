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

using BoSSS.Foundation;
using BoSSS.Solution;
using ilPSP.LinSolvers;
using NSE_SIMPLE.Multiphase;
using System.Collections.Generic;

namespace NSE_SIMPLE {

    /// <summary>
    /// Solver for level-set equation for multiphase flows.
    /// </summary>
    public class MultiphaseSolverLevelSet : SIMPLESolver {

        SIMPLEMatrixAssembly m_MatAsmblyLevelSet;
        SIMPLEMatrixAssembly m_MatAsmblyLevelSetApprox;
        ScalarFieldHistory<SinglePhaseField> m_LevelSet;
        double m_RelaxFactor;
        RelaxationTypes m_ModeRelaxLevelSet;

        // Time discretization
        BDFScheme m_BDF = null;

        MultiphaseSIMPLEControl m_multiphaseControl;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="solverConf"></param>
        /// <param name="sparseSolver"></param>        
        /// <param name="MatAsmblyScalar"></param>
        /// <param name="MatAsmblyScalarApprox"></param>
        /// <param name="Scalar"></param>        
        /// <param name="BDF"></param>
        public MultiphaseSolverLevelSet(SolverConfiguration solverConf, ISparseSolver sparseSolver,
            SIMPLEMatrixAssembly MatAsmblyScalar, SIMPLEMatrixAssembly MatAsmblyScalarApprox,
            ScalarFieldHistory<SinglePhaseField> Scalar, BDFScheme BDF)
            : base(solverConf, sparseSolver) {

            m_MatAsmblyLevelSet = MatAsmblyScalar;
            m_MatAsmblyLevelSetApprox = MatAsmblyScalarApprox;
            m_LevelSet = Scalar;

            m_solverConf = solverConf;
            m_multiphaseControl = solverConf.Control as MultiphaseSIMPLEControl;

            m_RelaxFactor = (1.0 - m_multiphaseControl.RelaxationFactorLevelSet) / m_multiphaseControl.RelaxationFactorLevelSet;
            m_ModeRelaxLevelSet = m_multiphaseControl.LevelSetRelaxationType;            

            m_BDF = BDF;
        }

        protected override MsrMatrix DefineMatrix(double dt) {
            MsrMatrix res = new MsrMatrix(m_MatAsmblyLevelSet.AssemblyMatrix);

            if ((m_ModeRelaxLevelSet == RelaxationTypes.Implicit) && (m_RelaxFactor != 0.0))
                res.Acc(m_RelaxFactor, m_MatAsmblyLevelSetApprox.AssemblyMatrix);

            if (m_BDF != null) {
                double LhsSummand = m_BDF.GetLhsSummand(dt, m_solverConf.BDFOrder);
                res.AccEyeSp(LhsSummand);
            }

            return res;
        }

        protected override IList<double> DefineRhs(double dt, int SpatialComponent) {
            double[] rhs = new double[m_MatAsmblyLevelSet.LocalLength];

            double[] LevelSetAffine = m_MatAsmblyLevelSet.AssemblyAffine;

            SinglePhaseField RhsSummand = new SinglePhaseField(m_LevelSet.Current.Basis);
            if (m_BDF != null)
                m_BDF.ComputeRhsSummand(dt, m_solverConf.BDFOrder, m_LevelSet, RhsSummand);

            double[] UnderRelaxation = new double[rhs.Length];
            if (((m_ModeRelaxLevelSet == RelaxationTypes.Implicit)) && (m_RelaxFactor != 0.0))
                m_MatAsmblyLevelSetApprox.AssemblyMatrix.SpMVpara(m_RelaxFactor, m_LevelSet.Current.CoordinateVector, 0.0, UnderRelaxation);

            for (int i = 0; i < rhs.Length; i++) {
                rhs[i] =
                    -LevelSetAffine[i]
                    - RhsSummand.CoordinateVector[i]
                    + UnderRelaxation[i];
            }

            return rhs;
        }
    }
}
