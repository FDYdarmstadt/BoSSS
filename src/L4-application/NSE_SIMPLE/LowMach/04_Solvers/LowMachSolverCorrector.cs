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
using BoSSS.Solution.NSECommon;

namespace NSE_SIMPLE {

    /// <summary>
    /// Corrector solver for Low-Mach flows.
    /// </summary>
    public class LowMachSolverCorrector : SIMPLESolver {

        SIMPLEMatrixAssembly MatAsmblyCorrector;

        SIMPLEOperator[] VelocityDivergence;
        VectorField<SinglePhaseField> Velocity_Intrmed;
        SinglePhaseField DivB4;

        BDFScheme BDF;
        ScalarFieldHistory<SinglePhaseField> Temperature;
        MaterialLaw EoS;

        double[] RHSManuDivKontiOperatorAffine = null;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="solverConf"></param>
        /// <param name="sparseSolver"></param>
        /// <param name="MatAsmblyCorrector"></param>
        /// <param name="VelocityDivergence"></param>
        /// <param name="Velocity_Intrmed"></param>
        /// <param name="DivB4"></param>
        /// <param name="BDF"></param>
        /// <param name="Temperature"></param>
        /// <param name="EoS"></param>
        public LowMachSolverCorrector(SolverConfiguration solverConf, ISparseSolver sparseSolver,
            SIMPLEMatrixAssembly MatAsmblyCorrector,
            SIMPLEOperator[] VelocityDivergence, VectorField<SinglePhaseField> Velocity_Intrmed, SinglePhaseField DivB4,
            BDFScheme BDF, ScalarFieldHistory<SinglePhaseField> Temperature, MaterialLaw EoS, double[] RHSManuDivKontiOperatorAffine = null)
            : base(solverConf, sparseSolver) {

            this.MatAsmblyCorrector = MatAsmblyCorrector;

            this.VelocityDivergence = VelocityDivergence;
            this.Velocity_Intrmed = Velocity_Intrmed;
            this.DivB4 = DivB4;

            this.BDF = BDF;
            this.Temperature = Temperature;
            this.EoS = EoS;

            this.RHSManuDivKontiOperatorAffine = RHSManuDivKontiOperatorAffine;
        }

        protected override MsrMatrix DefineMatrix(double dt) {
            MsrMatrix res = MatAsmblyCorrector.AssemblyMatrix;

            if (base.m_solverConf.Control.PressureReferencePoint != null)
                BoSSS.Solution.NSECommon.SolverUtils.SetRefPtPressure_Matrix(res, base.m_solverConf.PressureReferencePointIndex);

            return res;
        }

        protected override IList<double> DefineRhs(double dt, int SpatialComponent) {
            double[] Rhs = new double[DivB4.DOFLocal];

            
            // calculate mass defect
            DivB4.Clear();
            
            SolverUtils.CalculateMassDefect_Divergence(VelocityDivergence, Velocity_Intrmed, DivB4);
            Rhs.AccV(1.0, DivB4.CoordinateVector);
             
            // time derivative density
            if (base.m_solverConf.Control.Algorithm == SolutionAlgorithms.Unsteady_SIMPLE) {
                SinglePhaseField DensitySummand = new SinglePhaseField(DivB4.Basis);
                BDF.ComputeDensitySummand(dt, base.m_solverConf.BDFOrder, Temperature, EoS, DensitySummand);
                Rhs.AccV(1.0, DensitySummand.CoordinateVector);
            }

             
            
            // reference point pressure
            if (base.m_solverConf.Control.PressureReferencePoint != null)
                BoSSS.Solution.NSECommon.SolverUtils.SetRefPtPressure_Rhs(Rhs, base.m_solverConf.PressureReferencePointIndex, MatAsmblyCorrector.i0);

            return Rhs;
        }
    }
}
