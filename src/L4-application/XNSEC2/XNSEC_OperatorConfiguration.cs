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

using BoSSS.Application.XNSEC;
using BoSSS.Application.XNSFE_Solver;
using ilPSP;
using System;

namespace BoSSS.Solution.XNSECommon {

    public class XNSEC_OperatorConfiguration : XNSFE_OperatorConfiguration {

        public XNSEC_OperatorConfiguration(XNSEC_Control controlfile) : base(controlfile) {
            //m_ReactionRateConstants = controlfile.ReactionRateConstants;

            this.includeReactionTerms = controlfile.ChemicalReactionActive;
            this.manSolSource_OK = controlfile.ManufacturedSolutionSwitch;
            //this.physParams.IncludeConvection = controlfile.MomentumConvection_OK;
            //this.Transport = controlfile.MomentumConvection_OK;
            this.isSteady = (controlfile.TimesteppingMode == Control.AppControl._TimesteppingMode.Steady) ? true : false;

            this.timeDerivativeConti_OK = controlfile.timeDerivativeConti_OK;
            this.timeDerivativeEnergyp0_OK = controlfile.timeDerivativeEnergyp0_OK;

            this.Reynolds = controlfile.Reynolds;
            this.Prandtl = controlfile.Prandtl;
            this.Froude = controlfile.Froude;
            this.Lewis = controlfile.Lewis;
            this.NoOfChemicalSpecies = controlfile.NumberOfChemicalSpecies;
            this.gravityDirection = controlfile.GravityDirection;

            this.TemperatureEquationOK = controlfile.EnableTemperature;
            this.MassFractionEquationsOK = controlfile.EnableMassFractions;
            this.PlotAdditionalParameters = controlfile.PlotAdditionalParameters;

            //this.physParams = controlfile.PhysicalParametersCombustion;

            this.VariableReactionRateParameters = controlfile.VariableOneStepParameters;
            if (MassFractionEquationsOK == false && this.NoOfChemicalSpecies != 1)
                throw new Exception("Invalid configuration. ");
        }

        public bool VariableReactionRateParameters {
            get;
            set;
        }

        /// <summary>
        ///
        /// </summary>
        public bool timeDerivativeEnergyp0_OK {
            get;
            set;
        }

        public bool PlotAdditionalParameters {
            get;
            set;
        }

        public bool timeDerivativeConti_OK {
            get;
            set;
        }

        public bool includeReactionTerms {
            get;
            set;
        }

        /// <summary>
        /// If set to false, the temperature field will be set to 1.0 everywhere
        /// </summary>
        public bool TemperatureEquationOK {
            get;
            set;
        }

        /// <summary>
        /// If set to false, the mass fraction field 0 will be set equal to one.
        /// Note that
        /// </summary>
        public bool MassFractionEquationsOK {
            get;
            set;
        }

        public bool manSolSource_OK { get; set; }

        public bool isSteady { get; set; }

        public double[] Lewis {
            get;
            private set;
        }

        public double Reynolds {
            get;
            private set;
        }

        public double Prandtl {
            get;
            private set;
        }

        public double Froude {
            get;
            private set;
        }

        public int NoOfChemicalSpecies {
            get;
            private set;
        }

        public Vector gravityDirection {
            get;
            private set;
        }
    }
}