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
using ilPSP.LinSolvers;
using ilPSP.Utils;
using NSE_SIMPLE.Multiphase;
using System;
using System.Collections.Generic;

namespace NSE_SIMPLE {

    /// <summary>
    /// Corrector solver for multiphase flows.
    /// </summary>
    public class MultiphaseSolverCorrector : SIMPLESolver {

        SIMPLEOperator[] m_VelocityDivergence;
        SIMPLEMatrixAssembly m_MatAsmblyCorrector;

        VectorField<SinglePhaseField> m_Velocity_Intrmed;
        SinglePhaseField m_DivB4;

        /// <summary>
        /// Ctor without pressure stabilization.
        /// </summary>
        /// <param name="solverConfig"></param>
        /// <param name="_sparseSolver"></param>
        /// <param name="VelocityDivergence"></param>
        /// <param name="MatAsmblyCorrector"></param>                
        /// <param name="Velocity_Intrmed"></param>        
        /// <param name="DivB4"></param>        
        public MultiphaseSolverCorrector(SolverConfiguration solverConfig, ISparseSolver _sparseSolver,
            SIMPLEOperator[] VelocityDivergence, SIMPLEMatrixAssembly MatAsmblyCorrector,
            VectorField<SinglePhaseField> Velocity_Intrmed, SinglePhaseField DivB4)
            : base(solverConfig, _sparseSolver) {

            if ((VelocityDivergence.Length != solverConfig.SpatialDimension) || (Velocity_Intrmed.Dim != solverConfig.SpatialDimension))
                throw new ArgumentException("Mismatch of dimensions!");

            m_VelocityDivergence = VelocityDivergence;
            m_MatAsmblyCorrector = MatAsmblyCorrector;

            m_Velocity_Intrmed = Velocity_Intrmed;
            m_DivB4 = DivB4;
        }

        protected override MsrMatrix DefineMatrix(double dt) {
            MsrMatrix res = m_MatAsmblyCorrector.AssemblyMatrix;

            if (base.m_solverConf.Control.PressureReferencePoint != null)
                BoSSS.Solution.NSECommon.SolverUtils.SetRefPtPressure_Matrix(res, base.m_solverConf.PressureReferencePointIndex);

            return res;
        }

        protected override IList<double> DefineRhs(double dt, int SpatialComponent) {
            double[] Rhs = new double[m_DivB4.DOFLocal];

            // calculate mass defect
            m_DivB4.Clear();
            SolverUtils.CalculateMassDefect_Divergence(m_VelocityDivergence, m_Velocity_Intrmed, m_DivB4);
            Rhs.AccV(1.0, m_DivB4.CoordinateVector);

            // reference point pressure
            if (base.m_solverConf.Control.PressureReferencePoint != null)
                BoSSS.Solution.NSECommon.SolverUtils.SetRefPtPressure_Rhs(Rhs, base.m_solverConf.PressureReferencePointIndex, m_MatAsmblyCorrector.i0);

            return Rhs;
        }
    }
}
