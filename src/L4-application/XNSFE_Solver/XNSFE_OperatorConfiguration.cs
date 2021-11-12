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
using System.Threading.Tasks;

using ilPSP.Utils;

using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.EnergyCommon;
using ilPSP;
using BoSSS.Solution.NSECommon;
using BoSSS.Application.XNSE_Solver;

namespace BoSSS.Application.XNSFE_Solver {  
    public class XNSFE_OperatorConfiguration : XNSE_OperatorConfiguration, IXHeat_Configuration {


        public XNSFE_OperatorConfiguration(XNSFE_Control control)
            : base(control) {

            AgglomerationTreshold = control.AgglomerationThreshold;

            thermParams = control.ThermalParameters;

            solveEnergy = control.solveKineticEnergyEquation;

            HeatTransport = control.ThermalParameters.IncludeConvection;
            HeatSource = control.InitialValues_EvaluatorsVec.Keys.Any(name => name.StartsWith(VariableNames.HeatSource)) || control.FieldOptions.Keys.Where(k => k.Contains(VariableNames.HeatSource)).Any();
            solveHeat = control.solveCoupledHeatEquation;
            Evaporation = (control.ThermalParameters.hVap > 0.0);
            Buoyancy = control.ThermalParameters.alpha_A != 0.0 || control.ThermalParameters.alpha_B != 0.0;
            if(control.prescribedMassflux_Evaluator != null)
                prescribedMassflux = control.prescribedMassflux_Evaluator;
            if (control.prescribedMassflux != null)
                prescribedMassflux = control.prescribedMassflux.Evaluate;
            MatInt = !Evaporation;

            int nBlocks = 2;
            //if (solveEnergy) {
            //    CodBlocks = new bool[3];
            //    DomBlocks = new bool[3];
            //    nBlocks = 3;
            //}

            if (solveHeat) {
                this.conductMode = control.conductMode;
                CodBlocks = (this.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) ? new bool[nBlocks + 1] : new bool[nBlocks + 2];
                DomBlocks = (this.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP) ? new bool[nBlocks + 1] : new bool[nBlocks + 2];
            }

            CodBlocks.SetAll(true);
            DomBlocks.SetAll(true);
        }

        /// <summary>
        /// taken from <see cref="DoNotTouchParameters.CellAgglomerationThreshold"/>
        /// </summary>
        public double AgglomerationTreshold;

        /// <summary>
        /// 
        /// </summary>
        public ThermalParameters thermParams;

        /// <summary>
        /// include kinetic energy equation
        /// </summary>
        public bool solveEnergy;

        /// <summary>
        /// true if the heat equation is solved via the auxiliary heat flux formulation
        /// </summary>
        public ConductivityInSpeciesBulk.ConductivityMode conductMode;

        /// <summary>
        /// include heat equation
        /// </summary>
        public bool solveHeat;

        /// <summary>
        /// include transport operator
        /// </summary>
        public bool HeatTransport;

        /// <summary>
        /// include volumetric heat source
        /// </summary>
        public bool HeatSource;

        /// <summary>
        /// use upwind discretization
        /// </summary>
        public bool HeatUpwinding = false;

        /// <summary>
        /// include evaporation
        /// </summary>
        public bool Evaporation;

        /// <summary>
        /// include buoyancy, that is Boussinesq approximation
        /// </summary>
        public bool Buoyancy;

        /// <summary>
        /// 
        /// </summary>
        public Func<double[], double, double> prescribedMassflux;

        public ThermalParameters getThermParams {
            get { return thermParams; }
        }

        public ConductivityInSpeciesBulk.ConductivityMode getConductMode {
            get { return conductMode; }
        }

        public bool isHeatTransport {
            get { return HeatTransport; }
        }

        public bool isHeatSource {
            get { return HeatSource; }
        }

        public bool useUpwind {
            get { return HeatUpwinding;  }
        }

        public bool isEvaporation {
            get { return Evaporation; }
        }
        public bool isBuoyancy {
            get { return Buoyancy; }
        }
    }
}
