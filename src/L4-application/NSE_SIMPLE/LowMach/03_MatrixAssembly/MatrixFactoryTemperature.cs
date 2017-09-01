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
using BoSSS.Platform;
using NSE_SIMPLE.LowMach;

namespace NSE_SIMPLE {

    /// <summary>
    /// Collection of matrix assemblies for temperature equation.
    /// </summary>
    public class MatrixFactoryTemperature {

        /// <summary>
        /// Ctor.
        /// </summary>
        public MatrixFactoryTemperature(OperatorFactoryTemperature TemperatureOperators, int LocalNoOfCells, BlockDiagonalMatrix Rho,
            SolverConfiguration solverConf, BDFScheme BDF) {

            Temperature = new MatrixAssemblyTemperature(TemperatureOperators.TemperatureConvection, TemperatureOperators.HeatConduction);

            LowMachSIMPLEControl lowMachControl = solverConf.Control as LowMachSIMPLEControl;
            if (lowMachControl.RelaxationModeTemperature == RelaxationTypes.Implicit) {
                TemperatureApprox = new MatrixAssemblyApprox(
                    solverConf, LocalNoOfCells, Temperature, BDF, Rho, 2 * lowMachControl.PredictorApproximationUpdateCycle);
            }
        }

        /// <summary>
        /// Matrix for temperature equation, i.e. advection and heat conduction.
        /// </summary>
        public SIMPLEMatrixAssembly Temperature {
            get;
            private set;
        }

        /// <summary>
        /// Approximation of temperature matrix for under-relaxation.
        /// </summary>
        public SIMPLEMatrixAssembly TemperatureApprox {
            get;
            private set;
        }
    }
}
