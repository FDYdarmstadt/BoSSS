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

using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.LinSolvers;
using NSE_SIMPLE.BaseVariableDensity;
using System;

namespace NSE_SIMPLE.LowMach {

    /// <summary>
    /// Provides settings, which are special
    /// to the LowMach solver.
    /// </summary>
    public class LowMachSIMPLEControl : VariableDensitySIMPLEControl {

        public new MaterialLawLowMach EoS {
            get {
                return base.EoS as MaterialLawLowMach;
            }
            set {
                base.EoS = value;
            }
        }

        [ExclusiveLowerBound(0.0)]
        public double PenaltyHeatConductionScaling = 1.0;

        /// <summary>
        /// Under-relaxation factor for temperature.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        [InclusiveUpperBound(1.0)]
        public double RelexationFactorTemperature;

        /// <summary>
        /// Modus for under-relaxation of temperature.
        /// </summary>
        public RelaxationTypes RelaxationModeTemperature;

        /// <summary>
        /// Maximum number of iterations for flow solver within one segregated step.
        /// </summary>
        [InclusiveLowerBound(1)]
        public int MaxFlowSolverIterations;

        /// <summary>
        /// Maximum number of iterations for temperature solver within one segregated step.
        /// </summary>
        [InclusiveLowerBound(1)]
        public int MaxTemperatureSolverIterations;

        /// <summary>
        /// Convergence criterion for temperature.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double L2NormTemperatureResidual;

        /// <summary>
        /// Prandtl number.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double Prandtl;

        /// <summary>
        /// Heat capacity ratio
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double Gamma;

        /// <summary>
        /// Linear solver configuration temperature.
        /// </summary>
        [NotNull]
        public Func<ISparseSolver> TemperatureSolverFactory;

        /// <summary>
        /// Analytic solution temperature.
        /// </summary>
        public Func<double[], double> AnalyticTemperature = null;

        /// <summary>
        /// Modus for thermodynamic pressure.
        /// </summary>
        public ThermodynamicPressureMode ThermodynamicPressureMode;

        /// <summary>
        /// Has to be prescribed by the user for <see cref="ThermodynamicPressureMode.MassDetermined"/>.
        /// </summary>
        public double? InitialMass;

        /// <summary>
        /// EdgeTags for boundaries, where the Nusselt number
        /// should be calculated.
        /// </summary>
        public string[] EdgeTagsNusselt = null;

        /// <summary>
        /// 
        /// </summary>
        public override void Verify() {
            base.Verify();

            switch (ThermodynamicPressureMode) {
                case ThermodynamicPressureMode.Constant:
                    break;
                case ThermodynamicPressureMode.MassDetermined:
                    if (!InitialMass.HasValue) {
                        throw new Exception("InitialMass needs to be given in control-file, if Mode_ThermodynamicPressure is set to MassDetermined.");
                    }
                    break;
                default:
                    throw new NotImplementedException();
            }
        }
    }
}
