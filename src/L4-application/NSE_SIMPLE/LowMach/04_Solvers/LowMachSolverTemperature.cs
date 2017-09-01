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
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using NSE_SIMPLE.LowMach;

namespace NSE_SIMPLE {

    /// <summary>
    /// Solver for temperature equation for Low-Mach flows.
    /// </summary>
    public class LowMachSolverTemperature : SIMPLESolver {

        BlockDiagonalMatrix DensityMatrix;
        SIMPLEMatrixAssembly MatAsmblyTemperature;
        SIMPLEMatrixAssembly MatAsmblyTemperatureApprox;

        RelaxationTypes ModeRelaxTemperature;
        double RelaxFactor;
        ScalarFieldHistory<SinglePhaseField> Temperature;

        BDFScheme BDF;
        MaterialLaw EoS;

        double gamma;
        ScalarFieldHistory<SinglePhaseField> ThermodynamicPressure;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="solverConfig"></param>
        /// <param name="sparseSolver"></param>
        /// <param name="DensityMatrix"></param>
        /// <param name="MatAsmblyTemperature"></param>
        /// <param name="MatAsmblyTemperatureApprox"></param>
        /// <param name="Temperature"></param>
        /// <param name="BDF"></param>
        /// <param name="EoS"></param>
        /// <param name="ThermodynamicPressure"></param>
        public LowMachSolverTemperature(SolverConfiguration solverConfig, ISparseSolver sparseSolver,
            BlockDiagonalMatrix DensityMatrix, SIMPLEMatrixAssembly MatAsmblyTemperature, SIMPLEMatrixAssembly MatAsmblyTemperatureApprox,
            ScalarFieldHistory<SinglePhaseField> Temperature,
            BDFScheme BDF, MaterialLaw EoS,
            ScalarFieldHistory<SinglePhaseField> ThermodynamicPressure)
            : base(solverConfig, sparseSolver) {

            this.DensityMatrix = DensityMatrix;
            this.MatAsmblyTemperature = MatAsmblyTemperature;
            this.MatAsmblyTemperatureApprox = MatAsmblyTemperatureApprox;

            LowMachSIMPLEControl lowMachControl = solverConfig.Control as LowMachSIMPLEControl;

            this.ModeRelaxTemperature = lowMachControl.RelaxationModeTemperature;
            this.RelaxFactor = (1.0 - lowMachControl.RelexationFactorTemperature) / lowMachControl.RelexationFactorTemperature;
            this.Temperature = Temperature;

            this.BDF = BDF;
            this.EoS = EoS;

            this.gamma = lowMachControl.Gamma;
            this.ThermodynamicPressure = ThermodynamicPressure;
        }


        protected override MsrMatrix DefineMatrix(double dt) {
            MsrMatrix res = new MsrMatrix(MatAsmblyTemperature.AssemblyMatrix);

            if ((ModeRelaxTemperature == RelaxationTypes.Implicit) && (RelaxFactor != 0.0))
                res.Acc(RelaxFactor, MatAsmblyTemperatureApprox.AssemblyMatrix);

            if (BDF != null) {
                double LhsSummand = BDF.GetLhsSummand(dt, base.m_solverConf.BDFOrder);
                res.Acc(LhsSummand / gamma, DensityMatrix);
            }

            return res;
        }

        protected override IList<double> DefineRhs(double dt, int SpatialComponent) {
            double[] rhs = new double[MatAsmblyTemperature.LocalLength];

            double[] TemperatureAffine = MatAsmblyTemperature.AssemblyAffine;

            SinglePhaseField RhsSummand = new SinglePhaseField(Temperature.Current.Basis);
            if (base.m_solverConf.Control.Algorithm == SolutionAlgorithms.Unsteady_SIMPLE)
                BDF.ComputeRhsSummand(dt, base.m_solverConf.BDFOrder, Temperature, EoS, RhsSummand);

            double[] UnderRelaxation = new double[rhs.Length];
            if (((ModeRelaxTemperature == RelaxationTypes.Implicit)) && (RelaxFactor != 0.0))
                MatAsmblyTemperatureApprox.AssemblyMatrix.SpMVpara(RelaxFactor, Temperature.Current.CoordinateVector, 0.0, UnderRelaxation);

            for (int i = 0; i < rhs.Length; i++) {
                rhs[i] =
                    - TemperatureAffine[i]
                    - RhsSummand.CoordinateVector[i] / gamma                    
                    + UnderRelaxation[i];
            }

            return rhs;
        }
    }
}
