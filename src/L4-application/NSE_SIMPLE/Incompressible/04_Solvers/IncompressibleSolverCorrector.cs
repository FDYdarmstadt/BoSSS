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
using ilPSP.Utils;

namespace NSE_SIMPLE {

    /// <summary>
    /// Corrector solver for incompressible flows.
    /// </summary>
    public class SolverCorrector : SIMPLESolver {

        SolverConfiguration m_SolverConf;

        SIMPLEMatrixAssembly m_MatAsmblyCorrector;
        SIMPLEOperator[] m_VelocityDivergence;

        BDFScheme m_BDF = null;

        VectorField<SinglePhaseField> m_Velocity_Intrmed;
        SinglePhaseField m_DivB4;

        SIMPLEOperator m_PressureStabilization = null;
        SinglePhaseField m_Pressure;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="solverConf"></param>
        /// <param name="sparseSolver"></param>
        /// <param name="MatAsmblyCorrector"></param>
        /// <param name="VelocityDivergence"></param>
        /// <param name="BDF"></param>        
        /// <param name="Velocity_Intrmed"></param>
        /// <param name="DivB4"></param>
        /// <param name="PressureStabilization">Can be null</param>
        /// <param name="Pressure"></param>
        public SolverCorrector(SolverConfiguration solverConf, ISparseSolver sparseSolver,
            SIMPLEMatrixAssembly MatAsmblyCorrector, SIMPLEOperator[] VelocityDivergence, BDFScheme BDF,
            VectorField<SinglePhaseField> Velocity_Intrmed, SinglePhaseField DivB4,
            SIMPLEOperator PressureStabilization, SinglePhaseField Pressure)
            : base(solverConf, sparseSolver) {

            if ((VelocityDivergence.Length != solverConf.SpatialDimension) || (Velocity_Intrmed.Dim != solverConf.SpatialDimension))
                throw new ArgumentException("Mismatch of dimensions!");

            m_SolverConf = solverConf;

            m_MatAsmblyCorrector = MatAsmblyCorrector;
            m_VelocityDivergence = VelocityDivergence;

            m_BDF = BDF;

            m_Velocity_Intrmed = Velocity_Intrmed;
            m_DivB4 = DivB4;

            m_PressureStabilization = PressureStabilization;
            m_Pressure = Pressure;
        }

        /// <summary>
        /// Corrector matrix, cf. left-hand side of Eq. (17) in 
        /// B. Klein, F. Kummer, M. Keil, and M. Oberlack,
        /// An extension of the SIMPLE based discontinuous Galerkin solver to unsteady incompressible flows, J. Comput. Phys., 2013.
        /// </summary>
        /// <param name="dt"></param>
        /// <returns></returns>
        protected override MsrMatrix DefineMatrix(double dt) {
            MsrMatrix res = m_MatAsmblyCorrector.AssemblyMatrix;

            if (base.m_solverConf.Control.PressureReferencePoint != null)
                BoSSS.Solution.NSECommon.SolverUtils.SetRefPtPressure_Matrix(res, base.m_solverConf.PressureReferencePointIndex);

            return res;
        }

        /// <summary>
        /// Rhs corrector, cf. right-hand side of Eq. (17) in 
        /// B. Klein, F. Kummer, M. Keil, and M. Oberlack,
        /// An extension of the SIMPLE based discontinuous Galerkin solver to unsteady incompressible flows, J. Comput. Phys., 2013.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="SpatialComponent"></param>
        /// <returns></returns>
        protected override IList<double> DefineRhs(double dt, int SpatialComponent) {
            double[] Rhs = new double[m_DivB4.DOFLocal];

            // calculate mass defect
            m_DivB4.Clear();
            SolverUtils.CalculateMassDefect_Divergence(m_VelocityDivergence, m_Velocity_Intrmed, m_DivB4);
            Rhs.AccV(1.0, m_DivB4.CoordinateVector);

            // pressure stabilization
            if (base.m_solverConf.Control.PressureStabilizationScaling > 0.0) {
                double[] PressureStabi = new double[Rhs.Length];
                m_PressureStabilization.OperatorMatrix.SpMVpara(1.0, m_Pressure.CoordinateVector, 0.0, PressureStabi);
                Rhs.AccV(1.0, PressureStabi);
            }

            // reference point pressure
            if (base.m_solverConf.Control.PressureReferencePoint != null)
                BoSSS.Solution.NSECommon.SolverUtils.SetRefPtPressure_Rhs(Rhs, base.m_solverConf.PressureReferencePointIndex, m_MatAsmblyCorrector.i0);

            return Rhs;
        }
    }
}
