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
using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;

namespace NSE_SIMPLE {

    /// <summary>
    /// Collection of operators for level-set equation.
    /// </summary>
    public class OperatorFactoryLevelSet {

        /// <summary>
        /// Ctor.
        /// </summary>        
        /// <param name="ScalarMapping"></param>
        /// <param name="Velocity0"></param>
        /// <param name="Velocity0Mean"></param>        
        /// <param name="SolverConf"></param>
        public OperatorFactoryLevelSet(UnsetteledCoordinateMapping ScalarMapping,
            VectorField<SinglePhaseField> Velocity0, VectorField<SinglePhaseField> Velocity0Mean, SolverConfiguration SolverConf) {

            this.LevelSetAdvection = new LevelSetAdvectionOperator(ScalarMapping, Velocity0, Velocity0Mean, SolverConf);
        }        

        /// <summary>
        /// Advection operator for temperature
        /// </summary>
        public SIMPLEOperator LevelSetAdvection {
            get;
            private set;
        }
    }
}
